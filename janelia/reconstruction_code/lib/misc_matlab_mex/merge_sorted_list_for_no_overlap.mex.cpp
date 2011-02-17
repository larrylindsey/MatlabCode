// Given two sorted lists A, B with unique numbers compute a remapping
// of A such that it does not have any numbers in common with B.
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI.
//

#include <mex.h>
#include <math.h>
#include <algorithm>

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    mexPrintf("input params:\n");
    mexPrintf("\t1.	list a Px1 uint32\n");
    mexPrintf("\t2.	list b Qx1 uint32\n");
    mexPrintf("\tOutput:");
    mexPrintf("\t1. R X M*N matrix, where R is the number of samples in the input sampling mask.\n");
    mexPrintf("\t2. R X 2 matrix, coordinates of the samples on the image plane.\n");
    return;
  }
  if(nrhs!=2)
  {
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

  if(nlhs>1)
  {
    mexErrMsgTxt("Wrong number of outputs\n");
    return;
  }

	const mxArray * list_a_mx = prhs[0];
  const mxArray * list_b_mx = prhs[1];
  
  const int * size_a = mxGetDimensions(list_a_mx);
  int n_element_a = std::max(size_a[0], size_a[1]);
  const int * size_b = mxGetDimensions(list_b_mx);
  int n_element_b = std::max(size_b[0], size_b[1]);

  unsigned int * list_a = (unsigned int *) mxGetPr(list_a_mx);
  unsigned int * list_b = (unsigned int *) mxGetPr(list_b_mx);
  unsigned int max_element_a = 0;
  {
    int i;
    for(i=0; i<n_element_a; i++)
      max_element_a = std::max(max_element_a, list_a[i]);
  }
	plhs[0] = mxCreateNumericMatrix(max_element_a+1, 1, mxUINT32_CLASS,
                                  mxREAL);
  unsigned int * remap = (unsigned int*) mxGetPr(plhs[0]);
  
  // remap list A for non-overlapping merger with list B
  {
    int i, j=0;
    int current_offset = 1;
    for(i=0; i<n_element_a; i++){
      if(j>=n_element_b || i+current_offset<list_b[j]){
        remap[list_a[i]] = i + current_offset;
        continue;
      }
      while(j<n_element_b && i+current_offset==list_b[j]){
        current_offset++;
        j++;
      }
      remap[list_a[i]] = i + current_offset;      
    }
  }
  return;
}
