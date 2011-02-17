//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	08042008	init. code
//

#include <mex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ext/hash_map>
#include <utility>                       // for std::pair
#include <algorithm>                     // for std::for_each
#include <boost/graph/graph_traits.hpp>  // for boost::graph_traits

#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

namespace std{
  using namespace __gnu_cxx;
}

//                          !!IMPORTANT!!
// This must be the same as that used in the deformable mesh code
#define TRANSFORMATION_ID_OFFSET 10


struct equi{ bool operator()(const unsigned long int s1, const unsigned long int s2) const
    {return s1==s2;}};
typedef  std::hash_map<unsigned long int, unsigned int, std::hash<unsigned int>, equi> Hash_UInt64_UInt32;

using namespace boost;

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

	if(nrhs==0){
		if(nlhs==1){
			plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
			*mxGetPr(plhs[0]) = 1;
			return;
		}
		mexPrintf("Usage collect_seg_overlap_stats_area_paffine\n");
		mexPrintf("input params:\n");
		mexPrintf("\t1. superpixel label map 1 MxN matrix (uint32)\n");
		mexPrintf("\t2. superpixel label map 2 MxN matrix (uint32)\n");
		mexPrintf("\t3. scalar field 1 for computing histogram MxN matrix (uint8)\n");
		mexPrintf("\t4. scalar field 2 for computing histogram MxN matrix (uint8)\n");
    mexPrintf("\t5. piecewise mapping mask MxN uint16 matrix\n");
    mexPrintf("\t6. affine transformations 6xR double matrix\n");
		mexPrintf("output:\n");
		mexPrintf("\t1. label pair - area of overlap and sum of intensity Px(2+2) double.\n");
		return;
	}
	if(nrhs!=6){
		mexErrMsgTxt("Wrong number of inputs\n");
		return;
	}
	if(nlhs>1){
		mexErrMsgTxt("Wrong number of outputs\n");
		return;
	}

  mexPrintf("START: collect_seg_overlap_stats_interior_hist_pafffine\n");
	int numDim = mxGetNumberOfDimensions(prhs[0]);
	if(numDim!=2){
		mexErrMsgTxt("Wrong no. of dimensions for arg. 1\n");
		return;
	}

	const int * sizeImage;
	sizeImage = mxGetDimensions(prhs[0]);
	int width = sizeImage[1], height = sizeImage[0];
	int n_pixel = height*width;
	unsigned int * label_map_0 = (unsigned int *) mxGetPr(prhs[0]);
	unsigned int * label_map_1 = (unsigned int *) mxGetPr(prhs[1]);
	
	unsigned char * boundary_map_0 = (unsigned char *) mxGetPr(prhs[2]);
	unsigned char * boundary_map_1 = (unsigned char *) mxGetPr(prhs[3]);

  unsigned short int * map_mask = (unsigned short int *) mxGetPr(prhs[4]);
  double * transforms = mxGetPr(prhs[5]);
  
  // get sum and number of boundary values for each pair of adjacent segments
#ifdef __PROFILE__
  mexCallMATLAB(0, NULL, 0, NULL, "tic");
#endif

  Hash_UInt64_UInt32 label_pair_hist, overlap_interior_sum;
  {
    unsigned long int l, k;
    int transform_id;
    int x1, y1;
    for(int i=0, x=0; x<width; x++){
      for(int y=0; y<height; y++, i++){
        if(label_map_0[i]<=0)
          continue;
        if(map_mask[i]<TRANSFORMATION_ID_OFFSET)
          continue;
        transform_id = map_mask[i] - TRANSFORMATION_ID_OFFSET;
        x1 = (int) ((*(transforms + 6*transform_id))*x + 
                    (*(transforms + 6*transform_id + 2))*y +
                    (*(transforms + 6*transform_id + 4)));
        y1 = (int)((*(transforms + 6*transform_id + 1))*x + 
                   (*(transforms + 6*transform_id + 3))*y +
                   (*(transforms + 6*transform_id + 5)));
        if(x1<0 || x1>=width || y1<0 || y1>=height)
          continue;
        
        l = (unsigned long int) label_map_0[i];
        l = l << 32;
        k = (unsigned long int) label_map_1[x1*height+y1];
        l |= k;
        overlap_interior_sum[l] += (boundary_map_0[i]+boundary_map_1[x1*height+y1])/2;
        label_pair_hist[l]++;
      }
    }
  }
    
  {
    // output the histograms
    int n_label_pair = label_pair_hist.size();
    plhs[0] = mxCreateDoubleMatrix(n_label_pair, 2+2, mxREAL);
    double * overlap_area_out = mxGetPr(plhs[0]);
    int i=0;
    Hash_UInt64_UInt32::iterator lp_it;
    unsigned int s0, s1;
    unsigned long int lp, k, m;
    for(lp_it=label_pair_hist.begin(); lp_it!=label_pair_hist.end(); lp_it++){
      lp = (*lp_it).first;
      s0 = (unsigned int) (lp >> 32);
      s1 = (unsigned int) (lp & 0x00000000ffffffff);
      overlap_area_out[i] = (double) s0;
      overlap_area_out[i+n_label_pair] = (double) s1;
      overlap_area_out[i+2*n_label_pair] = (*lp_it).second;
      overlap_area_out[i+3*n_label_pair] = overlap_interior_sum[lp];
      i++;
    }
  }

  mexPrintf("STOP: collect_seg_overlap_stats_interior_hist_paffine\n");
  return;
}
