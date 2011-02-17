

#include <mex.h>
#include <math.h>
#include <stdio.h>

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
		mexPrintf("\t1.	NxN matrix of F\n");
		mexPrintf("\t2.	output file name\n");
    mexPrintf("output:\n");
    return;
  }
  if(nrhs!=2)
  {
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

	const mxArray * F_mx = prhs[0];
	const mxArray * output_filename_mx = prhs[1];
  
	const int * size_F = mxGetDimensions(F_mx);
	int n_node = size_F[0];
	double * F = mxGetPr(F_mx);
	
	char * output_filename = mxArrayToString(output_filename_mx);
  mexPrintf("output_filename: %s\n", output_filename);
  
  FILE * fout = fopen(output_filename, "wt");
  fprintf(fout, "min: ");
  
  for(int i=0; i<n_node; i++){
    for(int j=0; j<n_node; j++){
      fprintf(fout, "+ %g x_%d_%d ", F[i + j*n_node], i, j);
    }
  }
  fprintf(fout, ";\n");

  for(int i=0; i<n_node; i++){
    fprintf(fout, "x_%d_%d = 0;\n", i, i);
  }

  for(int i=0; i<n_node; i++){
    for(int j=0; j<n_node; j++){
      fprintf(fout, "0 <= x_%d_%d <= 1;\n", i, j);
      fprintf(fout, "x_%d_%d = x_%d_%d;\n", i, j, j, i);
    }
  }

  for(int i=0; i<n_node; i++){
    mexPrintf("i: %d\n", i);
    for(int j=0; j<n_node; j++){
      for(int k=0; k<n_node; k++){
        fprintf(fout, "x_%d_%d <= x_%d_%d + x_%d_%d;\n", i, k, i, j, j, k);
      }
    }
  }

  fclose(fout);
  
	mxFree(output_filename);
  
  return;
}
