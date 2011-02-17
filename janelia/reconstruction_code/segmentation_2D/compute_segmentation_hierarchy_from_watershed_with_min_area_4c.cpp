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

class Tree_Node{
	public:
		Tree_Node():	label(-1), parent(-1), child_0(-1), child_1(-1), f(-1), area(-1) {}
		int label;
    int label_new;
		int parent;
		int child_0, child_1;
		double f;
		int area;
};

// find most distant ancestor
int find_adam(Tree_Node * forest, int node_id){
	while(forest[node_id].parent!=-1)
		node_id = forest[node_id].parent;
	
	return node_id;
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
	if(numDim>2){
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
	Tree_Node * forest = (Tree_Node *) new Tree_Node[n_pixel*2];
	int p;
	double * ws_l = watershed_label_map;
	for(int i=0; i<n_pixel; i++, ws_l++){
		p=(int)*ws_l;
		if(p==0){
			n_boundary_elements++;
		}
		else{
			if(n_pixel_of_label[p]==0){ // first time this label is seen
				forest[n_label].label = p;
        forest[n_label].label_new = -1;
				forest[n_label].parent = -1;
				forest[n_label].child_0 = -1;
				forest[n_label].child_1 = -1;
				forest[n_label].f = -1;
				label_to_leaf[p]=n_label;
				n_label++;
			}
			n_pixel_of_label[p]++;
		}
	}	
	
	int n_node = n_label;
	for(int i=0; i<n_node; i++)
		forest[i].area = n_pixel_of_label[forest[i].label];
	
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
	double * f_map = watershed_field;
	Boundary * boundary_list = (Boundary *) new Boundary[n_boundary_elements*2];
	p_ws_l = padded_watershed_label_map+height+2+1;
	int j=0;
	int height_padded = height+2;
	for(int x=0; x<width; x++, p_ws_l+=2){
		for(int y=0; y<height; y++, f_map++, p_ws_l++){
			if((*p_ws_l)==0){
				boundary_list[j].f=*f_map; 	boundary_list[j].label_0=(int)*(p_ws_l-1); 						boundary_list[j].label_1=(int)*(p_ws_l+1);
        boundary_list[j].x = x; boundary_list[j].y = y;
				
        boundary_list[j+1].f=*f_map; 	boundary_list[j+1].label_0=(int)*(p_ws_l-height_padded); 		boundary_list[j+1].label_1=(int)*(p_ws_l+height_padded);
        boundary_list[j+1].x = x; boundary_list[j+1].y = y;

        j+=2;
			}
		}
	}
	
	// Get the list of boundary elements to be involved in the watershed
	std::sort(boundary_list, boundary_list+n_boundary_elements*2-1);
	
	Boundary_1 * boundary_list_sorted = (Boundary_1 *) new Boundary_1[n_boundary_elements*2];
	int length_boundary_list_sorted=0;
	Boundary * b_l = boundary_list;
	Boundary_1 *b_l_s = boundary_list_sorted;
	int prev_label_0=0, prev_label_1=0;
	for(int i=0; i<n_boundary_elements*2; i++, b_l++){
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
	
	std::sort(boundary_list_sorted, boundary_list_sorted+length_boundary_list_sorted-1);
	
	
	/////
	// merge the watershed tree
	/////
	int child_0_id, child_1_id;
	for(int i=0; i<length_boundary_list_sorted; i++){
		child_0_id = find_adam(forest, label_to_leaf[boundary_list_sorted[i].label_0]);
		child_1_id = find_adam(forest, label_to_leaf[boundary_list_sorted[i].label_1]);
		
		if(boundary_list_sorted[i].f > f_threshold){
			if( forest[child_0_id].area>=area_threshold && forest[child_1_id].area>=area_threshold )
				continue;
		}
		
		if(child_0_id==child_1_id)
			continue;
		
		int new_node_id=n_node;
		forest[new_node_id].label = n_label+n_node;
    forest[n_label].label_new = -1;
		forest[new_node_id].parent = -1;
		forest[new_node_id].child_0 = child_0_id;
		forest[new_node_id].child_1 = child_1_id;
		forest[new_node_id].f = boundary_list_sorted[i].f;
		forest[new_node_id].area = forest[child_0_id].area + forest[child_1_id].area;
		n_node++;
		
		forest[child_0_id].parent = new_node_id;
		forest[child_1_id].parent = new_node_id;
	}
	
	
	
	/////
	// Compute the new label map
	/////
	int * old_label_to_adam_node = (int *) new int[n_pixel];
	int * temp_node_list = (int *) new int[n_pixel];
	int tail;
	int n, adam;
  int n_new_label = 0;
	for(int i=0; i<n_node; i++){
		if(forest[i].parent!=-1)
			continue;
		adam = i;
		temp_node_list[0] = i;
    n_new_label++;
    forest[i].label_new = n_new_label;
		tail = 0;
		while(tail>=0){
			n = temp_node_list[tail];
			if(forest[n].child_0>=0){
				temp_node_list[tail+1] = forest[n].child_0;
				tail++;
				forest[n].child_0 = -2;
			}
			else{
				if(forest[n].child_1>=0){
					temp_node_list[tail+1] = forest[n].child_1;
					tail++;
					forest[n].child_1 = -2;
				}
				else{
					if(forest[n].child_0==-1 && forest[n].child_1==-1){
						old_label_to_adam_node[forest[n].label] = adam;
						tail--;
					}
					else{
						tail--;
					}
				}
			}
		}
	}
	
	plhs[0] = mxCreateDoubleMatrix(height, width, mxREAL);
	double * t_l = mxGetPr(plhs[0]);
	ws_l = watershed_label_map;
  int temp;
	for(int i=0; i<n_pixel; i++, ws_l++, t_l++){
    if(*ws_l==0)
      *t_l=0;
    else
      *t_l = forest[old_label_to_adam_node[(int)*ws_l]].label_new;
	}
	
	b_l = boundary_list;
  t_l = mxGetPr(plhs[0]);
	for(int i=0; i<n_boundary_elements*2; i++, b_l++){
		if(b_l->label_0<=0 || b_l->label_1<=0)
			continue;
    if(forest[old_label_to_adam_node[b_l->label_0]].label_new == forest[old_label_to_adam_node[b_l->label_1]].label_new)
      t_l[b_l->x*height + b_l->y] = forest[old_label_to_adam_node[b_l->label_0]].label_new;
  }  
	b_l = boundary_list;
	for(int i=0; i<n_boundary_elements*2; i++, b_l++){
		if(b_l->label_0<=0 || b_l->label_1<=0)
			continue;
    if(forest[old_label_to_adam_node[b_l->label_0]].label_new != forest[old_label_to_adam_node[b_l->label_1]].label_new)
      t_l[b_l->x*height + b_l->y] = 0;
  }  

  if(nlhs==2){
    const char * key_names[] = { "label", "parent", "child_0", "child_1", "f", "area"};
    plhs[1] = mxCreateStructMatrix(n_node, 1, 6, key_names);
    
    mxArray * output_forest = plhs[1];
    for(int i=0; i<n_node; i++){
      mxArray * v = mxCreateDoubleMatrix(1,1,mxREAL);
      *mxGetPr(v)=forest[i].label;
      mxSetFieldByNumber(output_forest, i, 0, v);
      
      v = mxCreateDoubleMatrix(1,1,mxREAL);
      *mxGetPr(v)=forest[i].parent!=-1?forest[i].parent+1:-1;
      mxSetFieldByNumber(output_forest, i, 1, v);
      
      v = mxCreateDoubleMatrix(1,1,mxREAL);
      *mxGetPr(v)=forest[i].child_0!=-1?forest[i].child_0+1:-1;
      mxSetFieldByNumber(output_forest, i, 2, v);
      
      v = mxCreateDoubleMatrix(1,1,mxREAL);
      *mxGetPr(v)=forest[i].child_1!=-1?forest[i].child_1+1:-1;
      mxSetFieldByNumber(output_forest, i, 3, v);
      
      v = mxCreateDoubleMatrix(1,1,mxREAL);
      *mxGetPr(v)=forest[i].f;
      mxSetFieldByNumber(output_forest, i, 4, v);
      
      v = mxCreateDoubleMatrix(1,1,mxREAL);
      *mxGetPr(v)=forest[i].area;
      mxSetFieldByNumber(output_forest, i, 5, v);
    }
  }
  
	delete [] n_pixel_of_label;
	delete [] label_to_leaf;
	delete [] forest;
	delete [] padded_watershed_label_map;
	delete [] boundary_list;
	delete [] boundary_list_sorted;
	delete [] old_label_to_adam_node;
	delete [] temp_node_list;
		
	return;
}
