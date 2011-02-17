// collect area of overlap between segments from two sections for linkage.
// A piecewise affine mapping is given from the first to the second segmentation map
// See align_stack_deformable_mesh_tile_pair_inter_plane.m
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// Matlab version called compute_segmentation_hierarchy_from_watershed_with_min_area.m
//
// v0	09282008	init. code
// v1 12122008  modified for non-prealigned segmentation maps. Use piecewise affine map
//                instead.
//

#include <mex.h>
#include <math.h>
#include <algorithm>
#include <ext/hash_map>

#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

#define MAX_3_3(p,H) MAX(MAX(MAX(MAX(*(p-1-H),*(p-1)),MAX(*(p-1+H),*(p+H))),MAX(MAX(*(p+1+H),*(p+1)),MAX(*(p+1-H),*(p-H)))),*(p))
#define MAX_2_2(p,H) MAX(MAX(*(p),*(p+1)),MAX(*(p+H),*(p+H+1))) 

//                          !!IMPORTANT!!
// This must be the same as that used in the deformable mesh code
#define TRANSFORMATION_ID_OFFSET 10

#define N_DIM 5

class Label_Pair{
	public:
		int label_0, label_1;
};
bool operator<(const Label_Pair& a, const Label_Pair& b){
	return (a.label_0<b.label_0)||((a.label_0==b.label_0)&&(a.label_1<b.label_1));
}

namespace std{
  using namespace __gnu_cxx;
}

struct equi{ bool operator()(const unsigned int s1, const unsigned int s2) const
    {return s1==s2;}};
typedef  std::hash_map<unsigned int, unsigned int,
                       std::hash<unsigned int>, equi> Hash_UInt32_UInt32;

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

	if(nrhs==0){
		mexPrintf("Usage features = collect_area_overlap_features_c(label_map_1, label_map_2);\n");
		mexPrintf("input params:\n");
		mexPrintf("\t1. label map 1 MxN matrix\n");
		mexPrintf("\t2. label map 2 MxN matrix\n");
    mexPrintf("\t3. piecewise mapping mask MxN uint16 matrix\n");
    mexPrintf("\t4. affine transformations 6xR double matrix\n");
		mexPrintf("output:\n");
		mexPrintf("\t1. features: [label 1, label 2, area of segment 1,\n");
		mexPrintf("\t\t area of segment 2, and area of overlap] Px5 matrix\n");
		return;
	}
	if(nrhs!=4){
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

  mexPrintf("START: collect_area_overlap_features_piecewise_affine_map_c\n");

  const mxArray * map_mask_mx = prhs[2];
  const mxArray * transforms_mx = prhs[3];
  
	const int * sizeImage;
	sizeImage = mxGetDimensions(prhs[0]);
	int width = sizeImage[1], height = sizeImage[0];
	int n_pixel = height*width;
	double * label_map_0 = mxGetPr(prhs[0]);
	double * label_map_1 = mxGetPr(prhs[1]);
  unsigned short int * map_mask = (unsigned short int *) mxGetPr(map_mask_mx);
  double * transforms = mxGetPr(transforms_mx);
  
	Hash_UInt32_UInt32 n_pixel_of_label_0;
  double * l;
  l = label_map_0;
	for(int i=0; i<n_pixel; i++, l++)
    if(*l>=0)
  		n_pixel_of_label_0[(unsigned int)*l]++;
	Hash_UInt32_UInt32 n_pixel_of_label_1;
  l = label_map_1;
	for(int i=0; i<n_pixel; i++, l++)
    if(*l>=0)
    	n_pixel_of_label_1[(int)*l]++;

  Label_Pair * label_pairs = (Label_Pair *) new Label_Pair[n_pixel];
  Label_Pair * l_p = label_pairs;
  double * l0 = label_map_0;
  unsigned short int * m_m = map_mask;
  double l1;
  int transform_id;
  int x1, y1;
  int n_pair = 0;
  for(int x=0; x<width; x++){
    for(int y=0; y<height; y++, l0++, m_m++){
      if(*l0<=0)
        continue;
      if(*m_m<TRANSFORMATION_ID_OFFSET)
        continue;
      transform_id = (*m_m) - TRANSFORMATION_ID_OFFSET;
      x1 = (int) ((*(transforms + 6*transform_id))*x + 
             (*(transforms + 6*transform_id + 2))*y +
             (*(transforms + 6*transform_id + 4)));
      y1 = (int)((*(transforms + 6*transform_id + 1))*x + 
             (*(transforms + 6*transform_id + 3))*y +
             (*(transforms + 6*transform_id + 5)));
      if(x1<0 || x1>=width || y1<0 || y1>=height)
        continue;
      
      l1 = (int)label_map_1[x1*height+y1];
      if(l1<=0)
        continue;
      l_p->label_0 = (int)*l0;
      l_p->label_1 = (int)l1;
      l_p++;
      n_pair++;
    }
  }
  mexPrintf("n_pair: %d\n", n_pair);
  
  if(n_pair==0){
    plhs[0] = mxCreateDoubleMatrix(N_DIM, 1, mxREAL);
    delete [] label_pairs;
    return;
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
  mexPrintf("n_unique_pair: %d\n", n_unique_pair);
  
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

  delete [] label_pairs;
  
  mexPrintf("STOP: collect_area_overlap_features_piecewise_affine_map_c\n");
	return;
}
