// Given two lists of segment-to-body mappings sorted according to
// body ids, generate pairs of segments that should be linked so that
// the segments form connected components corresponding to body ids.
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI.
//

#include <mex.h>
#include <math.h>
#include <algorithm>
#include <utility>
#include <vector>

typedef std::vector<std::pair<unsigned int, unsigned int> > Seg_Pair_Array;

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    mexPrintf("input params:\n");
    mexPrintf("\t1.	segment-to-body-map sorted by body ids Px2 uint32\n");
    mexPrintf("\t2.	segment-to-body-map sorted by body ids Qx2 uint32\n");
    mexPrintf("\tOutput:");
    mexPrintf("\t1. segment pairs to be linked Rx2 uint32\n");
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

	const mxArray * seg_2_body_1_mx = prhs[0];
  const mxArray * seg_2_body_2_mx = prhs[1];
  
  const int * size_1 = mxGetDimensions(seg_2_body_1_mx);
  int n_element_1 = std::max(size_1[0], size_1[1]);
  const int * size_2 = mxGetDimensions(seg_2_body_2_mx);
  int n_element_2 = std::max(size_2[0], size_2[1]);

  unsigned int * seg_2_body_1 = (unsigned int *) mxGetPr(seg_2_body_1_mx);
  unsigned int * seg_2_body_2 = (unsigned int *) mxGetPr(seg_2_body_2_mx);

  Seg_Pair_Array seg_pair;
  {
    std::pair<unsigned int, unsigned int> p;
    int i=0, j=0, h, k, i1, j1;
    unsigned int s_1, s_2, b_1, b_2, b_k, b_h;
    while(i<n_element_1 && j<n_element_2)
    {
      b_1 = seg_2_body_1[i+n_element_1];
      b_2 = seg_2_body_2[j+n_element_2];
      while(i<n_element_1 && b_1<b_2){
        i++;
        b_1 = seg_2_body_1[i+n_element_1];
      }
        while(j<n_element_2 && b_1>b_2){
          j++;
          b_2 = seg_2_body_2[j+n_element_2];
        }
        if(i>n_element_1 || j>=n_element_2)
          break;
        if(b_1!=b_2)
          continue;
        k=j;
        b_k = seg_2_body_2[k+n_element_2];
        while(k<n_element_2 && b_k==b_2){
          k++;
          b_k = seg_2_body_2[k+n_element_2];
        }
        h=i;
        b_h = seg_2_body_1[h+n_element_1];
        while(h<n_element_1 && b_h==b_1){
          h++;
          b_h = seg_2_body_1[h+n_element_1];
        }
        for(i1=i; i1<h; i1++)
          for(j1=j; j1<k; j1++){
            p.first = seg_2_body_1[i1];
            p.second = seg_2_body_2[j1];
            seg_pair.push_back(p);
          }
        i = h;
        j = k;
        }
    }
    
    plhs[0] = mxCreateNumericMatrix(seg_pair.size(), 2, mxUINT32_CLASS,
                                    mxREAL);
    {
      unsigned int n_pair = seg_pair.size();
      unsigned int * sp = (unsigned int*) mxGetPr(plhs[0]);
      int i=0;
      Seg_Pair_Array::iterator sp_it;
      for(sp_it=seg_pair.begin(); sp_it!=seg_pair.end(); sp_it++, i++){
        sp[i] = (*sp_it).first;
        sp[i+n_pair] = (*sp_it).second;
      }
    }
    return;
  }
