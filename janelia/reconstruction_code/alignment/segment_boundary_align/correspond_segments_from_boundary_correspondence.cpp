// Get correspondence between segments  of two segmentation maps given the
// correspondences between their boundary elements
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
//
// v0   08282008  init code
// v1   08312008  bug in label to vote alignment
//

#include <mex.h>
#include <math.h>
#include <algorithm>

#include <hash_functions.h>

typedef unsigned long int Label_Pair;
Label_Pair label_pair_2_id(unsigned int s1, unsigned int s2){
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


void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

	if(nrhs==0){
		if(nlhs==1){
			plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
			*mxGetPr(plhs[0]) = 1;
			return;
		}
		mexPrintf("Usage seg0_to_seg1_map =\n");
    mexPrintf("\tcorrespond_segments_from_boundary_correspondence(boundary0_to_boundary1_map, ...\n");
    mexPrintf("\t\tsegment_map_0, segment_map_1)\n");
		mexPrintf("input params:\n");
		mexPrintf("\t1. boundary0_to_boundary1_map: Map of boundary indexes of\n");
    mexPrintf("\t\tsegment_map_0 to segment_map_1 [2xN] matrix. Indexes as generated\n");
    mexPrintf("\t\tby sub2ind.\n");
		mexPrintf("\t2. segmentation_map_0 label map with 0's for boundaries - PxQ matrix.\n");
		mexPrintf("\t3. segmentation_map_1 PxQ matrix\n");
		mexPrintf("output:\n");
		mexPrintf("\t1. correspondence of labels of segment0 and segment1.\n\t\tEach column is [label_segment0, label_segment1, number of votes].\n");
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

	int numDim = mxGetNumberOfDimensions(prhs[1]);
	if(numDim!=2){
		mexErrMsgTxt("Wrong no. of dimensions for arg. 2\n");
		return;
	}

	double * b_2_b_map = mxGetPr(prhs[0]);
	double * segment_map_0 = mxGetPr(prhs[1]);
	double * segment_map_1 = mxGetPr(prhs[2]);
  
	const int * sizeImage;
	sizeImage = mxGetDimensions(prhs[1]);
	int width = sizeImage[1], height = sizeImage[0];
	int n_pixel = height*width;
	
  /* for debugging
  double a[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9};
  int b, c;
  get_two_labels(4, 3, 3, a, &b, &c);
  mexPrintf("%d %d\n", b, c);
  return;
  */

  Hash_UInt64_UInt32 label_pairs;
  
  const int * size_temp = mxGetDimensions(prhs[0]);
  int n_matched_boundary = size_temp[1];
  int index_0, index_1;
  int x_0, y_0, x_1, y_1;
  int n_segment_overlap_pair = 0;
  for(int i=0; i<2*n_matched_boundary; i+=2){
    index_0 = (int)b_2_b_map[i]-1; // 0 shifted
    index_1 = (int)b_2_b_map[i+1]-1;
    
    x_0 = index_0/height;
    y_0 = index_0%height;
    x_1 = index_1/height;
    y_1 = index_1%height;
    
    if(x_0>0 && x_0<width-1 && x_1>0 && x_1<width-1){
      if(segment_map_0[index_0-height]!=segment_map_0[index_0+height]){
        if(segment_map_1[index_1-height]!=segment_map_1[index_1+height]){
          label_pairs[label_pair_2_id((unsigned int)segment_map_0[index_0-height],
                  (unsigned int)segment_map_1[index_1-height])]++;
          label_pairs[label_pair_2_id((unsigned int)segment_map_0[index_0+height],
                  (unsigned int)segment_map_1[index_1+height])]++;
        }
      }
    }
    if(y_0>0 && y_0<height-1 && y_1>0 && y_1<height-1){
      if(segment_map_0[index_0-1]!=segment_map_0[index_0+1]){
        if(segment_map_1[index_1-1]!=segment_map_1[index_1+1]){
          label_pairs[label_pair_2_id((unsigned int)segment_map_0[index_0-1],
                  (unsigned int)segment_map_1[index_1-1])]++;
          label_pairs[label_pair_2_id((unsigned int)segment_map_0[index_0+1],
                  (unsigned int)segment_map_1[index_1+1])]++;
        }
      }
    }
  }
  
  #define N_OUT_DIM 3
	plhs[0] = mxCreateDoubleMatrix(N_OUT_DIM, label_pairs.size(), mxREAL);
  double * pairwise_votes = mxGetPr(plhs[0]);
  int i=0;
  for(Hash_UInt64_UInt32::iterator it = label_pairs.begin();
    it!=label_pairs.end(); it++, i+=N_OUT_DIM){
    std::pair<unsigned int, unsigned int> lp = id_2_label_pair((*it).first);
    pairwise_votes[i] = lp.first;
    pairwise_votes[i+1] = lp.second;
    pairwise_votes[i+2] = (*it).second;
  }
	
	return;
}

