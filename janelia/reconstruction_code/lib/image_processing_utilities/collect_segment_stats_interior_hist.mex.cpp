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

struct equi{ bool operator()(const unsigned int s1, const unsigned int s2) const
  {return s1==s2;}};
typedef  std::hash_map<unsigned int, unsigned int, std::hash<unsigned int>, equi> Hash_UInt32_UInt32;

// Hash function from pair of labels to unsigned int
struct equli{ bool operator()(const unsigned long int s1, const unsigned long int s2) const
  {return s1==s2;}};
typedef  std::hash_map<unsigned long int, unsigned int, std::hash<unsigned long int>, equli> Hash_UInt64_UInt32;

using namespace boost;

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    if(nlhs==1){
      plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
      *mxGetPr(plhs[0]) = 1;
      return;
    }
    mexPrintf("Usage collect_segment_pair_stats_boundary_hist\n");
    mexPrintf("input params:\n");
    mexPrintf("\t1. superpixel label map MxN matrix (uint32)\n");
    mexPrintf("\t2. scalar field on which segmentation was performed MxN matrix (uint8)\n");
    mexPrintf("\t3. bin size (uint8)\n");
    mexPrintf("output:\n");
    mexPrintf("\t1. label - histogram of interior values Pxn double.\n");
    return;
  }
  if(nrhs!=3){
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }
  if(nlhs>1){
    mexErrMsgTxt("Wrong number of outputs\n");
    return;
  }

  mexPrintf("START: collect_segment_stats_interior_hist\n");
  int numDim = mxGetNumberOfDimensions(prhs[0]);
  if(numDim!=2){
    mexErrMsgTxt("Wrong no. of dimensions for arg. 1\n");
    return;
  }

  const int * sizeImage;
  sizeImage = mxGetDimensions(prhs[0]);
  int width = sizeImage[1], height = sizeImage[0];
  int n_pixel = height*width;
  unsigned int * superpixel_label_map = (unsigned int *) mxGetPr(prhs[0]);
	
  unsigned char * boundary_map = (unsigned char *) mxGetPr(prhs[1]);

  unsigned int bin_size = (unsigned int) *((unsigned char *) mxGetPr(prhs[2]));

  // get sum and number of boundary values for each pair of adjacent segments
#ifdef __PROFILE__
  mexCallMATLAB(0, NULL, 0, NULL, "tic");
#endif

  unsigned int n_bin = (255/bin_size)+1;
  mexPrintf("n_bin: %d\n", n_bin);
  
  Hash_UInt64_UInt32 label_interior_hist;
  Hash_UInt32_UInt32 label_id_hist;
  {
    int i;
    unsigned long int l, k;
    for(i=0; i<n_pixel; i++){
      l = (unsigned long int) superpixel_label_map[i];
      l = l << 32;
      k = (unsigned long int) (boundary_map[i]/bin_size);
      l |= k;
      label_interior_hist[l]++;
      label_id_hist[superpixel_label_map[i]]++;
    }
  }
    
  {
    // output the histograms
    int n_label = label_id_hist.size();
    plhs[0] = mxCreateDoubleMatrix(n_label, n_bin+1, mxREAL);
    double * interior_hist_out = mxGetPr(plhs[0]);
    int i=0;
    Hash_UInt32_UInt32::iterator l_it;
    unsigned int l;
    unsigned long int k, m;
    for(l_it=label_id_hist.begin(); l_it!=label_id_hist.end(); l_it++){
      l = (*l_it).first;
      interior_hist_out[i] = l;
      for(int j=0; j<n_bin; j++){
        k = (unsigned long int) l;
        k = k << 32;
        m = (unsigned long int) j;
        k |= m;
        interior_hist_out[i+(1+j)*n_label] = (double)label_interior_hist[k];
      }      
      i++;
    }
  }

  mexPrintf("STOP: collect_segment_stats_interior_hist\n");
  return;
}
