// sum in all interval along row
// Shiv N. Vitaladevuni
//

#include <mex.h>
#include <math.h>

#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

/*
 * input params:
 *	1.	N x 2 matrix of N 2D vectors
 *	2.	1xp bins in 1st dimension
 *	3.	1xq bins in 2nd dimension
 *
 * output:
 *	Computes pxq 2D histogram of the vectors
 *
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
		mexPrintf("a = interval_sum(b)\n");
		mexPrintf("\tComputes sum in all intervals along rows.\n");
    mexPrintf("input params:\n");
		mexPrintf("\t1.	b: MxN matrix of M vectors\n");
    mexPrintf("output:\n");
		mexPrintf("\t1. a: Mx(N(N+1)/2) matrix\n");
    return;
  }
  if(nrhs!=1)
  {
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

	const mxArray * vectors_mx = prhs[0];
  int M = mxGetM(vectors_mx);
  int N = mxGetN(vectors_mx);

  // compute cumulative sum
  double * cumu_sum = (double *) new double[M*N];
  {
    double * vectors = mxGetPr(vectors_mx);
    int i, j, s;
    for(i=0; i<M; i++){
      s=0;
      for(j=0; j<N; j++){
        cumu_sum[i+j*M] = vectors[i+j*M] + s;
        s = cumu_sum[i+j*M];
      }
    }
  }
  
	plhs[0] = mxCreateDoubleMatrix(M, N*(N+1)/2, mxREAL);
  {
    double * interval_sum = mxGetPr(plhs[0]);
    int i, j, k, l;
    for(i=0; i<M; i++){
      l=0;
      for(j=0; j<N; j++){
        interval_sum[i + l*M] = cumu_sum[i + j*M];
        l++;
        for(k=j+1; k<N; k++){
          interval_sum[i + l*M] = cumu_sum[i + k*M]-cumu_sum[i + j*M];
          l++;
        }
      }
    }
  }

  delete [] cumu_sum;
  return;
}
