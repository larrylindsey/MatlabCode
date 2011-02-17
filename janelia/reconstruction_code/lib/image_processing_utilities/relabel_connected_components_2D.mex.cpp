// Relabel connected components of each segment with a distinct label
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
//

#include <mex.h>
#include <stdio.h>
#include <queue>

namespace std{
  using namespace __gnu_cxx;
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    if(nlhs==1){
      plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
      *mxGetPr(plhs[0]) = 1;
      return;
    }
    mexPrintf("Usage relabel_connected_components_2D\n");
    mexPrintf("input params:\n");
    mexPrintf("\t1. seg MxN uint32\n");
    mexPrintf("output:\n");
    mexPrintf("\t1. seg relabelled MxN uint32\n");
    return;
  }
  if(nrhs!=1){
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }
  if(nlhs>1){
    mexErrMsgTxt("Wrong number of outputs\n");
    return;
  }

  mexPrintf("START: relabel_connected_components_2D\n");
  const mxArray * seg_old_mx = prhs[0];
  const int M = mxGetM(seg_old_mx);
  const int N = mxGetN(seg_old_mx);
  int n_pixel = M*N;
  mexPrintf("M: %d, N: %d, n_pixel: %d\n", M, N, n_pixel);
  mexEvalString("drawnow;");
  unsigned int * seg_old_input = (unsigned int *) mxGetPr(seg_old_mx);

  unsigned int * seg_old = (unsigned int *) new unsigned int[n_pixel];
  {
    mexPrintf("copying seg map to temp. buffer\n");
    mexEvalString("drawnow;");
    int i;
    for(i=0; i<n_pixel; i++)
      seg_old[i] = seg_old_input[i];
  }

  plhs[0] = mxCreateNumericMatrix(M, N, mxUINT32_CLASS, mxREAL);
  unsigned int * seg_new = (unsigned int *) mxGetPr(plhs[0]);

  mexPrintf("relabelling\n");
  mexEvalString("drawnow;");
  int old_label_first_pos, pos;
  unsigned int old_label, new_label;
  std::queue<int> Q;
  old_label_first_pos = 0;
  new_label = 0;
  while(old_label_first_pos<n_pixel){
    old_label = seg_old[old_label_first_pos];
    if(old_label==0){
      // skip to next pixel
      old_label_first_pos++;
      continue;
    }

    // label assigned to all pixel in this connected component
    new_label++;

    // flood fill in this connected component
    Q.push(old_label_first_pos);
    seg_old[old_label_first_pos] = 0;
    while(!Q.empty()){
      pos = Q.front();
      Q.pop();
      seg_new[pos] = new_label;

      // enqueue this pixel's neighbors if have label = old_label
      int x = pos / M;
      int y = pos % M;

      if(x>0){ // east
	int pos_n = pos-M;
	if(seg_old[pos_n]==old_label){
	  seg_old[pos_n] = 0;
	  Q.push(pos_n);
	}
      }
      if(x<N-1){ // west
	int pos_n = pos+M;
	if(seg_old[pos_n]==old_label){
	  seg_old[pos_n] = 0;
	  Q.push(pos_n);
	}
      }
      if(y>0){ // north
	int pos_n = pos-1;
	if(seg_old[pos_n]==old_label){
	  seg_old[pos_n] = 0;
	  Q.push(pos_n);
	}
      }
      if(y<M-1){ // south
	int pos_n = pos+1;
	if(seg_old[pos_n]==old_label){
	  seg_old[pos_n] = 0;
	  Q.push(pos_n);
	}
      }
    }

    old_label_first_pos++;
  }

  delete [] seg_old;
  mexPrintf("STOP: relabel_connected_components_2D\n");
  return;
}
