// compute segmentation from a superpixel map (agglomerative
// clustering) Feature: between each segment, compute the mean
// boundary value and merge them if this value is below a specified
// threshold
// 
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
typedef  std::hash_map<unsigned int, unsigned int,
                       std::hash<unsigned int>, equi> Hash_UInt32_UInt32;

// Hash function from pair of labels to unsigned int
typedef unsigned long int Label_Pair;
struct equli{ bool operator()(const Label_Pair s1, const Label_Pair s2) const
    {return s1==s2;}};
typedef  std::hash_map<Label_Pair, unsigned int, std::hash<Label_Pair>, equli> Label_Pair_Hash;

Label_Pair label_pair_2_id(unsigned int s1, unsigned int s2){
  if(s1>s2){
    Label_Pair lp1;
    lp1 = s2;
    lp1 <<= 32;
    Label_Pair lp;
    lp = s1;
    lp |= lp1;
    return lp;
  }
  Label_Pair lp1;
  lp1 = s1;
  lp1 <<= 32;
  Label_Pair lp;
  lp = s2;
  lp |= lp1;
  return lp;
}
std::pair<unsigned int, unsigned int> id_2_label_pair(Label_Pair lp){
  std::pair<unsigned int, unsigned int> s;
  s.first = lp >> 32;
  s.second = lp & 0x00000000ffffffff;
  return s;
}

using namespace boost;

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

	if(nrhs==0){
		if(nlhs==1){
			plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
			*mxGetPr(plhs[0]) = 1;
			return;
		}
		mexPrintf("Usage old_label_to_new_label = segmentation_from_superpixels_mean_boundary_value_RAG(superpixel_label_map, boundary, f_threshold_seq, boundary_length_seq);\n");
		mexPrintf("input params:\n");
		mexPrintf("\t1. superpixel label map MxN matrix (uint32)\n");
		mexPrintf("\t2. scalar field on which segmentation was performed MxN matrix (uint8)\n");
		mexPrintf("output:\n");
		mexPrintf("\t1. label pair - sum and number of boundary values Px4 double.\n");
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

  mexPrintf("START: collect_segment_pair_stats_mean_boundary\n");
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

  // get sum and number of boundary values for each pair of adjacent segments
#ifdef __PROFILE__
  mexCallMATLAB(0, NULL, 0, NULL, "tic");
#endif
  Label_Pair_Hash sum_boundary_value, n_boundary_pixel;
  {
    int x, y;
    unsigned int * s;
    unsigned char * b;
    Label_Pair lp;
    for(x=1, s=superpixel_label_map+height, b=boundary_map+height;
        x<width-1;
        x++, s+=height, b+=height){
      if(*s==0){
        lp = label_pair_2_id(*(s-height), *(s+height));
        sum_boundary_value[lp] += *b;
        n_boundary_pixel[lp] ++;
      }
    }
    for(x=1, s=superpixel_label_map+height*2-1, b=boundary_map+height*2-1;
        x<width-1;
        x++, s+=height, b+=height){
      if(*s==0){
        lp = label_pair_2_id(*(s-height), *(s+height));
        sum_boundary_value[lp] += *b;
        n_boundary_pixel[lp] ++;
      }
    }
    for(y=1, s=superpixel_label_map+1, b=boundary_map+1;
        y<height-1;
        y++, s++, b++){
      if(*s==0){
        lp = label_pair_2_id(*(s-1), *(s+1));
        sum_boundary_value[lp] += *b;
        n_boundary_pixel[lp] ++;
      }
    }
    for(y=1, s=superpixel_label_map+height*(width-1)+1, b=boundary_map+height*(width-1)+1;
        y<height-1;
        y++, s++, b++){
      if(*s==0){
        lp = label_pair_2_id(*(s-1), *(s+1));
        sum_boundary_value[lp] += *b;
        n_boundary_pixel[lp] ++;
      }
    }
    for(x=1, s=superpixel_label_map+height, b=boundary_map+height; // skip first column
        x<width-1;
        x++){
      s++; // skip first element in column
      b++;
      for(y=1;
          y<height-1;
          y++, s++, b++){
        if(*s==0){
          lp = label_pair_2_id(*(s-height), *(s+height));
          sum_boundary_value[lp] += *b;
          n_boundary_pixel[lp] ++;
        }
        if(*s==0){
          lp = label_pair_2_id(*(s-1), *(s+1));
          sum_boundary_value[lp] += *b;
          n_boundary_pixel[lp] ++;
        }
        if(*s==0){
          lp = label_pair_2_id(*(s-height-1), *(s+height+1));
          sum_boundary_value[lp] += *b;
          n_boundary_pixel[lp] ++;
        }
        if(*s==0){
          lp = label_pair_2_id(*(s-height+1), *(s+height-1));
          sum_boundary_value[lp] += *b;
          n_boundary_pixel[lp] ++;
        }
      }
      s++; // skip last element in column.
      b++;
    }
    
#ifdef __DEBUG__
    {
      Label_Pair_Hash::iterator lp_it;
      unsigned int s0, s1;
      Label_Pair lp;
      for(lp_it=sum_boundary_value.begin();
          lp_it!=sum_boundary_value.end(); lp_it++){
        lp = (*lp_it).first;
        tie(s0,s1) = id_2_label_pair(lp);
        mexPrintf("label0: %u, label1: %u, sum_boundary_value: %u, n_boundary_pixel: %u\n",
          s0, s1, (*lp_it).second, n_boundary_pixel[lp]);
      }
    }
#endif
    
  }

  {
    // output the mean boundary values
    int n_pair = sum_boundary_value.size();
    plhs[0] = mxCreateDoubleMatrix(n_pair, 4, mxREAL);
    double * mean_boundary = mxGetPr(plhs[0]);
    int i=0;
    Label_Pair_Hash::iterator lp_it;
    Label_Pair lp;
    unsigned int s0, s1;
    for(lp_it=sum_boundary_value.begin();
        lp_it!=sum_boundary_value.end(); lp_it++){
      lp = (*lp_it).first;
      tie(s0,s1) = id_2_label_pair(lp);
      if(s0==0 || s1==0)
        continue;
      mean_boundary[i] = s0;
      mean_boundary[i+n_pair] = s1;
      mean_boundary[i+2*n_pair] = (double)sum_boundary_value[lp];
      mean_boundary[i+3*n_pair] = (double)n_boundary_pixel[lp];
      i++;
    }
  }

  mexPrintf("STOP: collect_segment_pair_stats_mean_boundary\n");
  return;
}
