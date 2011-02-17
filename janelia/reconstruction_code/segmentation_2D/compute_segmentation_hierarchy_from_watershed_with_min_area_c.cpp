// compute watershed hierarchy such that each segment has atleast minimum area
// An alternative to Ladder.m designed by Yuriy Mischenko ~ Dec. 2007
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// Matlab version called compute_segmentation_hierarchy_from_watershed_with_min_area.m
//
// v0	01042008	init. code
//

#include <mex.h>
#include <math.h>
#include <algorithm>
#include <ext/hash_map>

#include <merge_sets_h.h>

namespace std{
  using namespace __gnu_cxx;
}

struct eqi{ bool operator()(const int s1, const int s2) const
    {return s1==s2;}};

#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

#define MAX_3_3(p,H) MAX(MAX(MAX(MAX(*(p-1-H),*(p-1)),MAX(*(p-1+H),*(p+H))),MAX(MAX(*(p+1+H),*(p+1)),MAX(*(p+1-H),*(p-H)))),*(p))
#define MAX_2_2(p,H) MAX(MAX(*(p),*(p+1)),MAX(*(p+H),*(p+H+1))) 

#define NUM_ORIENT_BINS		180

class Boundary{
	public:
		int label_0, label_1;
    int x,y;
    double f;
};
bool operator<(const Boundary& a, const Boundary& b){
	return (a.label_0<b.label_0)||((a.label_0==b.label_0)&&(a.label_1<b.label_1))||((a.label_0==b.label_0)&&(a.label_1==b.label_1)&&(a.f<b.f));
}

class Boundary_1{
	public:
		int label_0, label_1;
		double f;
};
bool operator<(const Boundary_1& a, const Boundary_1& b){
	return a.f<b.f;
}

/*
 * input params:
 * 	1. watershed label map MxN matrix
 *	2. scalar field on which watershed was performed MxN matrix
 *	3. threshold on the scalar field below which all watershed segments must be merged
 *	4. threshold on the area - minimum number pixels in a segment for it to exist on its own. If violated, the segment is merged.
 *
 * output:
 *	1. merged watershed label map
 *	2. watershed merged hierarchy - a forest
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

	if(nrhs==0){
		if(nlhs==1){
			plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
			*mxGetPr(plhs[0]) = 1;
			return;
		}
		mexPrintf("Usage [watershed_merged_hierarchy, label_map] = compute_segmentation_hierarchy_from_watershed_with_min_area_c(watershed_img_c, boundary, f_threshold, area_threshold);\n");
		mexPrintf("input params:\n");
		mexPrintf("\t1. watershed label map MxN matrix\n");
		mexPrintf("\t2. scalar field on which watershed was performed MxN matrix\n");
		mexPrintf("\t3. threshold on the scalar field below which all watershed segments must be merged\n");
		mexPrintf("\t4. threshold on the area - minimum number pixels in a segment for it to exist on its own. If violated, the segment is merged.\n");
		mexPrintf("output:\n");
		mexPrintf("\t1. merged watershed label map\n");
		mexPrintf("\t2. watershed merged hierarchy - a forest\n");
		return;
	}
	if(nrhs!=4){
		mexErrMsgTxt("Wrong number of inputs\n");
		return;
	}


	int numDim = mxGetNumberOfDimensions(prhs[0]);
	if(numDim!=2){
		mexErrMsgTxt("Wrong no. of dimensions for arg. 1\n");
		return;
	}

	const int * sizeImage;
	sizeImage = mxGetDimensions(prhs[0]);
	int width = sizeImage[1], height = sizeImage[0];
	int n_pixel = height*width;
	double * watershed_label_map = mxGetPr(prhs[0]);
	
	double * watershed_field = mxGetPr(prhs[1]);
	double f_threshold = *(mxGetPr(prhs[2]));
	int area_threshold = (int) *(mxGetPr(prhs[3]));
	
	int * n_pixel_of_label = (int *) new int[n_pixel];
	int * l=n_pixel_of_label;
	for(int i=0; i<n_pixel; i++, l++)
		*l=0;
	
	int n_boundary_elements=0;
	int n_label=0;
	int * label_to_leaf = (int *) new int[n_pixel];
	int p;
	double * ws_l = watershed_label_map;
	for(int i=0; i<n_pixel; i++, ws_l++){
		p=(int)*ws_l;
		if(p==0)
			n_boundary_elements++;
		else
			n_pixel_of_label[p]++;
	}	
	
	/////
	// Construct a list of boundary elements
	/////
	// Padding to take care of image edges
	double * padded_watershed_label_map = (double *) new double[(height+2)*(width+2)];
	double * p_ws_l = padded_watershed_label_map;
	for(int i=0; i<(height+2)*(width+2); i++, p_ws_l++)
		*p_ws_l=0;
		
	p_ws_l = padded_watershed_label_map+height+2+1;
	ws_l = watershed_label_map;
	for(int x=0; x<width; x++, p_ws_l+=2){
		for(int y=0; y<height; y++, ws_l++,p_ws_l++){
			*p_ws_l=*ws_l;
		}
	}
	
	// Collect elements
  Merge_Sets_H<int, std::hash<int>, eqi>
    superpixel_sets(NULL);
	double * f_map = watershed_field;
	Boundary * boundary_list = (Boundary *) new Boundary[n_boundary_elements*4];
	p_ws_l = padded_watershed_label_map+height+2+1;
	int j=0;
	int height_padded = height+2;
  int label_0, label_1;
	for(int x=0; x<width; x++, p_ws_l+=2){
		for(int y=0; y<height; y++, f_map++, p_ws_l++){
			if((*p_ws_l)==0){
        label_0=(int)*(p_ws_l-1);
        label_1=(int)*(p_ws_l+1);
        if(label_0!=0 && label_1!=0){
          if(label_0<label_1){
            boundary_list[j].label_0=label_0;
            boundary_list[j].label_1=label_1;
          }
          else{
            boundary_list[j].label_0=label_1;
            boundary_list[j].label_1=label_0;
          }
          boundary_list[j].f=*f_map;
          boundary_list[j].x = x;
          boundary_list[j].y = y;
          j++;
        }
        
        label_0=(int)*(p_ws_l-height_padded);
        label_1=(int)*(p_ws_l+height_padded);
        if(label_0!=0 && label_1!=0){
          if(label_0<label_1){
            boundary_list[j].label_0=label_0;
            boundary_list[j].label_1=label_1;
          }
          else{
            boundary_list[j].label_0=label_1;
            boundary_list[j].label_1=label_0;
          }
          boundary_list[j].f=*f_map;
          boundary_list[j].x = x;
          boundary_list[j].y = y;
          j++;
        }
        /*
        label_0=(int)*(p_ws_l-height_padded-1);
        label_1=(int)*(p_ws_l+height_padded+1);
        if(label_0!=0 && label_1!=0){
          if(label_0<label_1){
            boundary_list[j].label_0=label_0;
            boundary_list[j].label_1=label_1;
          }
          else{
            boundary_list[j].label_0=label_1;
            boundary_list[j].label_1=label_0;
          }
          boundary_list[j].f=*f_map;
          boundary_list[j].x = x;
          boundary_list[j].y = y;
          j++;
        }
        
        label_0=(int)*(p_ws_l-height_padded+1);
        label_1=(int)*(p_ws_l+height_padded-1);
        if(label_0!=0 && label_1!=0){
          if(label_0<label_1){
            boundary_list[j].label_0=label_0;
            boundary_list[j].label_1=label_1;
          }
          else{
            boundary_list[j].label_0=label_1;
            boundary_list[j].label_1=label_0;
          }
          boundary_list[j].f=*f_map;
          boundary_list[j].x = x;
          boundary_list[j].y = y;
          j++;
        }
        */
			}
      else{
        superpixel_sets.add_new_set_inc(*p_ws_l);
      }
		}
	}
  n_boundary_elements = j;
  
  mexPrintf("n_boundary_elements: %d\n", n_boundary_elements);
  
  if(n_boundary_elements==0){
    plhs[0] = mxCreateDoubleMatrix(height, width, mxREAL);
    double * t_l = mxGetPr(plhs[0]);
    ws_l = watershed_label_map;
    int temp;
    for(int i=0; i<n_pixel; i++, ws_l++, t_l++)
      *t_l=*ws_l;
    delete [] n_pixel_of_label;
    delete [] label_to_leaf;
    delete [] padded_watershed_label_map;
    delete [] boundary_list;
    return;
  }  
	
	// Get the list of boundary elements to be involved in the watershed
	std::sort(boundary_list, boundary_list+n_boundary_elements-1);
	
	Boundary_1 * boundary_list_sorted = (Boundary_1 *) new Boundary_1[n_boundary_elements];
	int length_boundary_list_sorted=0;
	Boundary * b_l = boundary_list;
	Boundary_1 *b_l_s = boundary_list_sorted;
	int prev_label_0=0, prev_label_1=0;
	for(int i=0; i<n_boundary_elements; i++, b_l++){
		if(b_l->label_0==0 || b_l->label_1==0)
			continue;
		
		if(b_l->label_1!=prev_label_1 || b_l->label_0!=prev_label_0){
			b_l_s->label_0 = b_l->label_0;
			b_l_s->label_1 = b_l->label_1;
			b_l_s->f = b_l->f;
			b_l_s++;
			length_boundary_list_sorted++;
			
			prev_label_0=b_l->label_0;
			prev_label_1=b_l->label_1;
		}
	}
	
  mexPrintf("length_boundary_list_sorted: %d\n", length_boundary_list_sorted);
  
  if(length_boundary_list_sorted==0){
    plhs[0] = mxCreateDoubleMatrix(height, width, mxREAL);
    double * t_l = mxGetPr(plhs[0]);
    ws_l = watershed_label_map;
    int temp;
    for(int i=0; i<n_pixel; i++, ws_l++, t_l++)
      *t_l=*ws_l;
    delete [] n_pixel_of_label;
    delete [] label_to_leaf;
    delete [] padded_watershed_label_map;
    delete [] boundary_list;
    delete [] boundary_list_sorted;
    return;
  }  
	
	std::sort(boundary_list_sorted, boundary_list_sorted+length_boundary_list_sorted-1);
	mexEvalString("drawnow;");
  
	/////
	// merge the watershed tree
	/////
	int child_0_id, child_1_id;
	for(int i=0; i<length_boundary_list_sorted; i++){
    /*
    if(i%10000==0){
      mexPrintf("i: %d\n", i);
      mexEvalString("drawnow;");
    }
    */
		child_0_id = superpixel_sets.get_adam(boundary_list_sorted[i].label_0);
		child_1_id = superpixel_sets.get_adam(boundary_list_sorted[i].label_1);
		
		if(child_0_id==child_1_id)
			continue;
		
		if(boundary_list_sorted[i].f > f_threshold){
			if(n_pixel_of_label[child_0_id]>=area_threshold &&
         n_pixel_of_label[child_1_id]>=area_threshold )
				continue;
		}
		
    superpixel_sets.merge(child_0_id, child_1_id);
    n_pixel_of_label[child_0_id] += n_pixel_of_label[child_1_id];
    n_pixel_of_label[child_1_id] = n_pixel_of_label[child_0_id];
	}
	
	plhs[0] = mxCreateDoubleMatrix(height, width, mxREAL);
	double * t_l = mxGetPr(plhs[0]);
	ws_l = watershed_label_map;
  int temp;
	for(int i=0; i<n_pixel; i++, ws_l++, t_l++){
    if(*ws_l==0)
      *t_l=0;
    else
      *t_l = superpixel_sets.get_adam(*ws_l);
	}

  /*
  b_l = boundary_list;
  t_l = mxGetPr(plhs[0]);
	for(int i=0; i<n_boundary_elements; i++, b_l++){
		if(b_l->label_0<=0 || b_l->label_1<=0)
			continue;
    if(superpixel_sets.get_adam(b_l->label_0) ==
       superpixel_sets.get_adam(b_l->label_1))
      t_l[b_l->x*height + b_l->y] = superpixel_sets.get_adam(b_l->label_1);
  }  
  
	b_l = boundary_list;
	for(int i=0; i<n_boundary_elements; i++, b_l++){
		if(b_l->label_0<=0 || b_l->label_1<=0)
			continue;
    if(superpixel_sets.get_adam(b_l->label_0) !=
       superpixel_sets.get_adam(b_l->label_1))
      t_l[b_l->x*height + b_l->y] = 0;
  }
  */ 
  
	delete [] n_pixel_of_label;
	delete [] label_to_leaf;
	delete [] padded_watershed_label_map;
	delete [] boundary_list;
	delete [] boundary_list_sorted;
		
	return;
}
