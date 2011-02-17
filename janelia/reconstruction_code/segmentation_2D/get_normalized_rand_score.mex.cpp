// compute normalized rand score for evaluating 2D segmentation.
// Rand index is normalized by the area of the two segments
// Shiv N. Vitaladevuni
//

#include <mex.h>
#include <math.h>
#include <hash_functions.h>

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    mexPrintf("input params:\n");
    mexPrintf("\t1. seg 1 uint32\n");
    mexPrintf("\t2. seg 2 uint32\n");
    mexPrintf("\t3. skip size uint32\n");
    mexPrintf("output:\n");
    mexPrintf("\t [rand_score_false_merge, rand_score_false_split]\n");
    return;
  }
  if(nrhs!=3){
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

  mexPrintf("START: get_normalized_rand_score\n");

  const mxArray * seg_1_mx = prhs[0];
  const mxArray * seg_2_mx = prhs[1];
  const mxArray * skip_size_mx = prhs[2];

  const int M = mxGetM(seg_1_mx);
  const int N = mxGetN(seg_1_mx);
  mexPrintf("M: %d, N: %d\n", M, N);
  const unsigned int skip_size = * (unsigned int *) mxGetPr(skip_size_mx);
  unsigned int * seg_1 = (unsigned int *) mxGetPr(seg_1_mx);
  unsigned int * seg_2 = (unsigned int *) mxGetPr(seg_2_mx);
  mexEvalString("drawnow;");

  mexPrintf("computing area of segments\n");
  Hash_Int32_Int32 area_1, area_2;
  {
    int i;
    for(i=0; i<M*N; i++){
      area_1[seg_1[i]]++;
      area_2[seg_2[i]]++;
    }
  }

  int xa, ya, xb, yb;
  double rand_score_false_merge=0, rand_score_false_split=0, rs_fm, rs_fs;
  bool p1, p2;
  int ida, idb;
  unsigned int a1, a2;
  double na;
  for(xa=0; xa<N; xa+=skip_size){
    for(ya=0; ya<M; ya+=skip_size){
      ida = xa*M + ya;

      if(seg_1[ida]==0 || seg_2[ida]==0)
	continue;

      rs_fm = 0;
      rs_fs = 0;
      for(xb=0; xb<N; xb+=skip_size){
	for(yb=0; yb<M; yb+=skip_size){
	  idb = xb*M + yb;

	  if(seg_1[idb]==0 || seg_2[idb]==0)
	    continue;

	  p1 = seg_1[ida] == seg_1[idb];
	  p2 = seg_2[ida] == seg_2[idb];
#ifdef __DEBUG__
	  mexPrintf("%u %u %u %u %d\n", seg_1[ida], seg_1[idb], seg_2[ida],
		    seg_2[idb], (p1 && !p2) || (!p1 && p2));
	  mexEvalString("drawnow;");
#endif
	  rs_fm += (!p1 && p2)/(double)area_1[seg_1[idb]];
	  rs_fs += (p1 && !p2)/(double)area_1[seg_1[idb]];
	}
      }
      rs_fm /= (double) area_1[seg_1[ida]];
      rs_fs /= (double) area_1[seg_1[ida]];

      rand_score_false_merge += rs_fm;
      rand_score_false_split += rs_fs;
    }
  }  
  plhs[0] = mxCreateDoubleMatrix(1, 2, mxREAL);
  mxGetPr(plhs[0])[0] = rand_score_false_merge;
  mxGetPr(plhs[0])[1] = rand_score_false_split;

  mexPrintf("STOP: get_normalized_rand_score\n");
  mexEvalString("drawnow;");
  return;
}
