// collect area of overlap between segments from two sections for linkage.
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// Matlab version called compute_segmentation_hierarchy_from_watershed_with_min_area.m
//
// v0	09282008	init. code
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

class Label_Pair{
	public:
		int label_0, label_1;
};
bool operator<(const Label_Pair& a, const Label_Pair& b){
	return (a.label_0<b.label_0)||((a.label_0==b.label_0)&&(a.label_1<b.label_1));
}


/*
 * input params:
 * 	1. label map 1 MxN matrix
 *	2. label map 2 MxN matrix
 *
 * output:
 *	1. features: [label 1, label 2, area of segment 1, 
 *                area of segment 2, and area of overlap] Px5 matrix 
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

	if(nrhs==0){
		mexPrintf("Usage features = collect_area_overlap_features_c(label_map_1, label_map_2);\n");
		mexPrintf("input params:\n");
		mexPrintf("\t1. label map 1 MxN matrix\n");
		mexPrintf("\t2. label map 2 MxN matrix\n");
		mexPrintf("output:\n");
		mexPrintf("\t1. features: [label 1, label 2, area of segment 1,\n");
		mexPrintf("\t\t area of segment 2, and area of overlap] Px5 matrix\n");
		return;
	}
	if(nrhs!=2){
		mexErrMsgTxt("Wrong number of inputs\n");
		return;
	}
	if(nlhs>1){
		mexErrMsgTxt("Wrong number of outputs\n");
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
	double * label_map_0 = mxGetPr(prhs[0]);
	double * label_map_1 = mxGetPr(prhs[1]);
	
	int * n_pixel_of_label_0 = (int *) new int[n_pixel];
	int * n = n_pixel_of_label_0;
	for(int i=0; i<n_pixel; i++, n++)
		*n=0;
  double * l;
  l = label_map_0;
	for(int i=0; i<n_pixel; i++, l++)
    if(*l>=0)
  		n_pixel_of_label_0[(int)*l]++;
	int * n_pixel_of_label_1 = (int *) new int[n_pixel];
	n = n_pixel_of_label_1;
	for(int i=0; i<n_pixel; i++, n++)
		*n=0;
  l = label_map_1;
	for(int i=0; i<n_pixel; i++, l++)
    if(*l>=0)
    	n_pixel_of_label_1[(int)*l]++;

  Label_Pair * label_pairs = (Label_Pair *) new Label_Pair[n_pixel];
  Label_Pair * l_p = label_pairs;
  double * l0 = label_map_0, * l1 = label_map_1;
  int n_pair = 0;
	for(int i=0; i<n_pixel; i++, l0++, l1++){
    if(*l0<=0 || *l1<=0)
      continue;
    
    l_p->label_0 = (int)*l0;
    l_p->label_1 = (int)*l1;
    l_p++;
    n_pair++;
  }

	// sort the label pairs
	std::sort(label_pairs, label_pairs+n_pair-1);
	
  // count number of overlapping segment pairs
	Label_Pair * label_pairs_sorted = label_pairs;
	int n_unique_pair=0;
	l_p = label_pairs_sorted;
	int prev_label_0=0, prev_label_1=0;
	for(int i=0; i<n_pair; i++, l_p++){
		if(l_p->label_0<=0 || l_p->label_1<=0)
			continue;
		
		if(l_p->label_1!=prev_label_1 || l_p->label_0!=prev_label_0){
      n_unique_pair++;
			prev_label_0=l_p->label_0;
			prev_label_1=l_p->label_1;
		}
	}
  
  #define N_DIM 5
  plhs[0] = mxCreateDoubleMatrix(N_DIM, n_unique_pair, mxREAL);
  
  // compute features
	l_p = label_pairs_sorted;
	prev_label_0=0;
  prev_label_1=0;
  int count = 0;
  double * feature = mxGetPr(plhs[0]);
	for(int i=0; i<n_pair; i++, l_p++){
		if(l_p->label_0<=0 || l_p->label_1<=0)
			continue;
		
		if(l_p->label_1!=prev_label_1 || l_p->label_0!=prev_label_0){
      if(count!=0){
        *(feature+4)  = count;
        feature += N_DIM;
      }
      
      *feature      = l_p->label_0;
      *(feature+1)  = l_p->label_1;
      *(feature+2)  = n_pixel_of_label_0[l_p->label_0];
      *(feature+3)  = n_pixel_of_label_1[l_p->label_1];
      count = 1;
      
			prev_label_0=l_p->label_0;
			prev_label_1=l_p->label_1;
		}
    else
      count ++;
	}
  *(feature+4)  = count;

  delete [] n_pixel_of_label_0;
  delete [] n_pixel_of_label_1;
  delete [] label_pairs;
  
	return;
}
