// compute histogram of neighboring intensities
// for mitochondria detection
// Shiv N. Vitaladevuni
//

#include <mex.h>
#include <math.h>

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    mexPrintf("input params:\n");
    mexPrintf("\t1. seg 1 uint32\n");
    mexPrintf("\t2. seg 2 uint32\n");
    mexPrintf("\t3. skip size uint32\n");
    mexPrintf("output:\n");
    mexPrintf("\t[rand_score_false_merge, rand_score_false_split]\n");
    return;
  }
  if(nrhs!=3)
    {
      mexErrMsgTxt("Wrong number of inputs\n");
      return;
    }

  mexPrintf("start: get_rand_score2\n");

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

  int xa, ya, xb, yb;
  unsigned int rand_score_false_merge=0, rand_score_false_split=0;
  bool p1, p2;
  unsigned int * s1a, * s1b, * s2a, * s2b;
  for(xa=0; xa<N; xa+=skip_size){
    s1a = seg_1 + xa*M;
    s2a = seg_2 + xa*M;
    for(ya=0; ya<M; ya+=skip_size, s1a+=skip_size, s2a+=skip_size){
      if(*s1a==0 || *s2a==0)
	continue;
      for(xb=0; xb<N; xb+=skip_size){
	s1b = seg_1 + xb*M;
	s2b = seg_2 + xb*M;
	for(yb=0; yb<M; yb+=skip_size, s1b+=skip_size, s2b+=skip_size){
	  if(*s2a==0 || *s2b==0)
	    continue;

	  p1 = *s1a == *s1b;
	  p2 = *s2a == *s2b;
#ifdef __DEBUG__
	  mexPrintf("%u %u %u %u %d\n", *s1a, *s1b, *s2a, *s2b,
		    (p1 && !p2) || (!p1 && p2));
	  mexEvalString("drawnow;");
#endif
	  rand_score_false_merge += !p1 && p2;
	  rand_score_false_split += p1 && !p2;
	}
      }
    }
  }  
  plhs[0] = mxCreateDoubleMatrix(1, 2, mxREAL);
  mxGetPr(plhs[0])[0] = (double) rand_score_false_merge;
  mxGetPr(plhs[0])[1] = (double) rand_score_false_split;

  mexPrintf("stop: get_rand_score2\n");
  mexEvalString("drawnow;");
  return;
}
