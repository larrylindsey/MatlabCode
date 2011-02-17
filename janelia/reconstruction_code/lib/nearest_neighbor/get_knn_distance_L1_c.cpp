// compute nearest neighbor L1 distances of a set of vectors to  another set
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	10242008	init. code
//

#include <mex.h>
#include <math.h>
#include <algorithm>

#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

#define MAX_3_3(p,H) MAX(MAX(MAX(MAX(*(p-1-H),*(p-1)),MAX(*(p-1+H),*(p+H))),MAX(MAX(*(p+1+H),*(p+1)),MAX(*(p+1-H),*(p-H)))),*(p))
#define MAX_2_2(p,H) MAX(MAX(*(p),*(p+1)),MAX(*(p+H),*(p+H+1))) 

#define NUM_ORIENT_BINS		180
#define INF 1000000

/*
 * input params:
 * 	1. feature vector set A: PxM matrix of M feature vectors
 *	2. feature vector set B: PxN matrix of N feature vectors
 *
 * output:
 *	1. distances 1xM: distance to nearest vector in set B
 *  2. min. id. 1xM: id of the nearest vector in set B
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

	if(nrhs==0){
		if(nlhs==1){
			plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
			*mxGetPr(plhs[0]) = 1;
			return;
		}
		mexPrintf("Usage [distance, min_id] = get_knn_distance_L1_c(feature_set_A, feature_set_B, k);\n");
		mexPrintf("input params:\n");
		mexPrintf("\t1. feature vector set A: PxM matrix of M feature vectors\n");
		mexPrintf("\t2. feature vector set B: PxN matrix of N feature vectors\n");
		mexPrintf("\t3. k: natural number\n");
		mexPrintf("output:\n");
		mexPrintf("\t1. distances Mx1: distance to k nearest vectors in set B\n");
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

	const int * size_0 = mxGetDimensions(prhs[0]);
	int n_vect_0 = size_0[1], n_dim_0 = size_0[0];
	const int * size_1 = mxGetDimensions(prhs[1]);
	int n_vect_1 = size_1[1], n_dim_1 = size_1[0];
	if(n_dim_0!=n_dim_1){
		mexErrMsgTxt("Wrong number of dimensions for the two sets of vectors\n");
		return;
	}
  int n_dim = n_dim_0;
  
	double * vects_0 = mxGetPr(prhs[0]);
	double * vects_1 = mxGetPr(prhs[1]);
  int k_value = (int) * mxGetPr(prhs[2]);
	k_value = MIN(k_value, n_vect_1);
  
	plhs[0] = mxCreateDoubleMatrix(n_vect_0, 1, mxREAL);
  double * knn_distance = mxGetPr(plhs[0]);
  double * k_d = knn_distance;
  
  double * distances = (double *) new double[n_vect_1];
  double * d_p;
  
  double d;
  double * v0, * v1;
  double * vects_0_p;
  vects_0_p = vects_0;
  for(int i=0; i<n_vect_0; i++, vects_0_p+=n_dim, k_d++){
    v1 = vects_1;
    d_p = distances;
    for(int j=0; j<n_vect_1; j++, d_p++){
      v0 = vects_0_p;
      d=0;
      for(int k=0; k<n_dim; k++, v0++, v1++){
        d += fabs((*v0) - (*v1));
      }

      *d_p = d;
    }
    std::sort(distances, distances + n_vect_1);
    *k_d = 0;
    d_p = distances;
    for(int j=0; j<k_value; j++, d_p++)
      *k_d += *d_p;
    *k_d /= (double) k_value;
  }
  
  delete [] distances;
	return;
}
