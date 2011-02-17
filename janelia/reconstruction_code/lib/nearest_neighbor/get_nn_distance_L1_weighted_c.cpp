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
 *	3. feature weights PxN matrix of N weight vectors
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
		mexPrintf("Usage [distance, min_id] = get_fn_distance_L1_weighted_c(\n");
    mexPrintf("\t\tfeature_set_A, feature_set_B, weights);\n");
		mexPrintf("input params:\n");
		mexPrintf("\t1. feature vector set A: PxM matrix of M feature vectors\n");
		mexPrintf("\t2. feature vector set B: PxN matrix of N feature vectors\n");
		mexPrintf("\t3. feature weights PxN matrix of N weight vectors\n");
		mexPrintf("output:\n");
		mexPrintf("\t1. distances Mx1: distance to nearest vector in set B\n");
		mexPrintf("\t2. min. id. Mx1: id of the nearest vector in set B\n");
		return;
	}
	if(nrhs!=3){
		mexErrMsgTxt("Wrong number of inputs\n");
		return;
	}
	if(nlhs>2){
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
  
	const int * size_w = mxGetDimensions(prhs[2]);
	int n_vect_w = size_w[1], n_dim_w = size_w[0];
	if(n_dim_w!=n_dim || n_vect_w!=n_vect_1){
		mexErrMsgTxt("Dimensions of set B and weights don't match\n");
		return;
	}
  
	double * vects_0 = mxGetPr(prhs[0]);
	double * vects_1 = mxGetPr(prhs[1]);
	double * weights = mxGetPr(prhs[2]);
	
	plhs[0] = mxCreateDoubleMatrix(n_vect_0, 1, mxREAL);
  double * min_distance = mxGetPr(plhs[0]);
  double * m_d = min_distance;
  for(int i=0; i<n_vect_0; i++, m_d++)
    *m_d = INF;
  
  double * m_i, * min_id;
  if(nlhs==2){
  	plhs[1] = mxCreateDoubleMatrix(n_vect_0, 1, mxREAL);
    min_id = mxGetPr(plhs[1]);
  }
  
  if(nlhs==2){
    double d;
    double * v0, * v1, * w;
    double * vects_1_p, * weights_p;
    vects_1_p = vects_1;
    weights_p = weights;
    for(int i=0; i<n_vect_1; i++, vects_1_p+=n_dim, weights_p+=n_dim){
      v0 = vects_0;
      m_d = min_distance;
      m_i = min_id;
      for(int j=0; j<n_vect_0; j++, m_d++, m_i++){
        v1 = vects_1_p;
        w = weights_p;
        d=0;
        for(int k=0; k<n_dim; k++, v0++, v1++, w++){
          d += (*w) * fabs((*v0) - (*v1));
        }
        if(d<(*m_d)){
          *m_d = d;
          *m_i = i+1; // in MATLAB's indexing
        }
      }
    }
  }else{
    double d;
    double * v0, * v1, * w;
    double * vects_1_p, * weights_p;
    vects_1_p = vects_1;
    weights_p = weights;
    for(int i=0; i<n_vect_1; i++, vects_1_p+=n_dim, weights_p+=n_dim){
      v0 = vects_0;
      m_d = min_distance;
      for(int j=0; j<n_vect_0; j++, m_d++){
        v1 = vects_1_p;
        w = weights_p;
        d=0;
        for(int k=0; k<n_dim; k++, v0++, v1++, w++){
          d += (*w) * fabs((*v0) - (*v1));
        }
        if(d<(*m_d)){
          *m_d = d;
        }
      }
    }
  }
	return;
}
