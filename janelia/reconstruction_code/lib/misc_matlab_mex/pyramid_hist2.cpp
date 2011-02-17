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
		mexPrintf("\t1.	MxN 2D matrix\n");
		mexPrintf("\t2.	scale factor s - kth level is M/(s^k) x N/(s^k)\n");
    mexPrintf("output:\n");
		mexPrintf("\t1D vector: concatenation of multiresolution pyramid of 2D matrix\n");
    return;
  }
  if(nrhs!=2){
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

	const mxArray * matrix_mx = prhs[0];
	
	const int scale_factor = (int) * mxGetPr(prhs[1]);
	if(scale_factor<2)
	{
		mexErrMsgTxt("Scale factor has to be greater than 1\n");
		return;
	}
  
	const int * size_matrix = mxGetDimensions(matrix_mx);
	int M = size_matrix[0];
	int N = size_matrix[1];
	double * matrix = mxGetPr(matrix_mx);
	
	int M1 = M+1;
	int N1 = N+1;
	double * integral_matrix = (double *) new double[M1*N1];
	double sum=0;
	for(int j=0; j<M1; j++)
		integral_matrix[j]=0;
	for(int i=0; i<N1; i++)
		integral_matrix[i*M1]=0;
	for(int j=0; j<M; j++){
		sum += matrix[j];
		integral_matrix[M1+j+1] = sum;
	}
	for(int i=1; i<N; i++){
		sum = matrix[i*M];
		integral_matrix[(i+1)*M1+1] = sum + integral_matrix[i*M1+1];
		for(int j=1; j<M; j++){
			sum += matrix[i*M+j];
			integral_matrix[(i+1)*M1+j+1] = sum + integral_matrix[i*M1+j+1];
		}
	}
	
	int skip = scale_factor;
	double * buffer = (double *) new double[M*N*3];
	int buffer_id=0;
	while(skip<=M && skip<=N){
		for(int i=skip; i<N1; i+=skip){
			for(int j=skip; j<M1; j+=skip){
				buffer[buffer_id] = integral_matrix[i*M1+j] + integral_matrix[(i-skip)*M1+j-skip] - integral_matrix[(i-skip)*M1+j] - integral_matrix[i*M1+j-skip];
				buffer_id++;
			}
			if(M%skip!=0){
				buffer[buffer_id] = integral_matrix[i*M1+M] + integral_matrix[(i-skip)*M1+M-skip] - integral_matrix[(i-skip)*M1+M] - integral_matrix[i*M1+M-skip];
				buffer_id++;
			}
		}
		if(N%skip!=0){
			for(int j=skip; j<M1; j+=skip){
				buffer[buffer_id] = integral_matrix[N*M1+j] + integral_matrix[(N-skip)*M1+j-skip] - integral_matrix[(N-skip)*M1+j] - integral_matrix[N*M1+j-skip];
				buffer_id++;
			}
			if(M%skip!=0){
				buffer[buffer_id] = integral_matrix[N*M1+M] + integral_matrix[(N-skip)*M1+M-skip] - integral_matrix[(N-skip)*M1+M] - integral_matrix[N*M1+M-skip];
				buffer_id++;
			}
		}
		
		skip *= scale_factor;
	}
	/*
	plhs[0] = mxCreateDoubleMatrix(M1, N1, mxREAL);
	double * pyramid_out = mxGetPr(plhs[0]);
	for(int i=0; i<M1*N1; i++)
		pyramid_out[i] = integral_matrix[i];
	*/
  
	plhs[0] = mxCreateDoubleMatrix(1, buffer_id, mxREAL);
	double * pyramid_out = mxGetPr(plhs[0]);
	for(int i=0; i<buffer_id; i++)
		pyramid_out[i] = buffer[i];
	
	return;
}
