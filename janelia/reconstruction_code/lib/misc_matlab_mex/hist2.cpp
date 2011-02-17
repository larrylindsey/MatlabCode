// compute histogram of neighboring intensities
// for mitochondria detection
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
    mexPrintf("input params:\n");
		mexPrintf("\t1.	Nx2 matrix of N 2D vectors\n");
		mexPrintf("\t2.	1xp bins in 1st dimension\n");
		mexPrintf("\t3.	1xq bins in 2nd dimension\n");
    mexPrintf("output:\n");
		mexPrintf("\tComputes pxq 2D histogram of the vectors\n");
    return;
  }
  if(nrhs!=3)
  {
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

	const mxArray * vectors_mx = prhs[0];
	const mxArray * bins_1_mx = prhs[1];
	const mxArray * bins_2_mx = prhs[2];
  
	const int * size_vectors = mxGetDimensions(vectors_mx);
	int n_vector = size_vectors[0];
	double * vectors = mxGetPr(vectors_mx);
	
	const int * size_bins_1 = mxGetDimensions(bins_1_mx);
	const int n_bin_1 = size_bins_1[0]*size_bins_1[1];
	double * bins_1 = mxGetPr(bins_1_mx);
	double * bin_boundaries_1 = (double *) new double[n_bin_1-1];
	for(int i=0; i<n_bin_1-1; i++){
		bin_boundaries_1[i] = (bins_1[i]+bins_1[i+1])/2.0;
	}
	
	const int * size_bins_2 = mxGetDimensions(bins_2_mx);
	const int n_bin_2 = size_bins_2[0]*size_bins_2[1];
	double * bins_2 = mxGetPr(bins_2_mx);
	double * bin_boundaries_2 = (double *) new double[n_bin_2-1];
	for(int i=0; i<n_bin_2-1; i++){
		bin_boundaries_2[i] = (bins_2[i]+bins_2[i+1])/2.0;
	}
	
	const int n_bin = n_bin_1 * n_bin_2;
	
	plhs[0] = mxCreateDoubleMatrix(n_bin_1, n_bin_2, mxREAL);
	double * hist = mxGetPr(plhs[0]);
	for(int i=0; i<n_bin; i++)
		hist[i]=0;
	
	int j, bin_id_1, bin_id_2;
	double * b_b_p;
	double * v_p_1 = vectors, * v_p_2 = vectors+n_vector;
	for(int i=0; i<n_vector; i++, v_p_1++, v_p_2++){
		// first dimension
		b_b_p = bin_boundaries_1;
		for(j=0; j<n_bin_1-1; j++, b_b_p++){
			if(*v_p_1<*b_b_p)
				break;
		}
		bin_id_1 = j;
	
		// second dimension
		b_b_p = bin_boundaries_2;
		for(j=0; j<n_bin_2-1; j++, b_b_p++){
			if(*v_p_2<*b_b_p)
				break;
		}
		bin_id_2 = j;
	
		(*(hist + bin_id_2*n_bin_1 + bin_id_1))++;
	}
			 
	delete [] bin_boundaries_1;
	delete [] bin_boundaries_2;
  
  return;
}
