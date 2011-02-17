// Compute all pairwise distances (Euclidean) among a set of vectors
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	12292008	init. code
//

#include <mex.h>
#include <math.h>
#include <algorithm>

#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

/*
 * input params:
 * 	1. feature vector set A: NxP matrix of N feature vectors
 *
 * output:
 *	1. distances 1 x (N*(N-1)/2): distance to nearest vector in set A
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

	if(nrhs==0){
		if(nlhs==1){
			plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
			*mxGetPr(plhs[0]) = 1;
			return;
		}
		mexPrintf("Usage [distance, min_id] = get_nn_distance_L1_c(feature_set_A, feature_set_B);\n");
		mexPrintf("input params:\n");
		mexPrintf("\t1. feature vector set A: NxP matrix of N feature vectors\n");
		mexPrintf("output:\n");
		mexPrintf("\t1. distances 1 x (N*(N-1)/2): distance to nearest vector in set A\n");
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

	const int * size_0 = mxGetDimensions(prhs[0]);
	int n_vect = size_0[0], n_dim = size_0[1];
  
	double * vects = mxGetPr(prhs[0]);
	
	plhs[0] = mxCreateDoubleMatrix(n_vect*(n_vect-1)/2, 1, mxREAL);
  double * distances = mxGetPr(plhs[0]);
  double * d_p;
  double d;
  double * v0, * v1;
  for(int k=0; k<n_dim; k++){
    v0 = vects + k*n_vect;
    d_p = distances;
    for(int i=0; i<n_vect-1; i++, v0++){
      v1 = v0 + 1;
      for(int j=i+1; j<n_vect; j++, v1++, d_p++){
        d = (*v0) - (*v1);
        *d_p += d*d;
      }
    }
  }
  
  d_p = distances;
  for(int i=0; i<n_vect*(n_vect-1)/2; i++, d_p++)
    *d_p = sqrt(*d_p);
  
	return;
}
