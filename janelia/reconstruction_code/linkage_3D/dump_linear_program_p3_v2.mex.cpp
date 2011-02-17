#include <mex.h>
#include <math.h>
#include <stdio.h>

#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    mexPrintf("input params:\n");
		mexPrintf("\t1.	NxN matrix of F\n");
		mexPrintf("\t2.	NxN adjacency matrix\n");
		mexPrintf("\t3.	output file name\n");
    return;
  }
  if(nrhs!=3)
  {
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

	const mxArray * F_mx = prhs[0];
	const mxArray * adjacency_mx = prhs[1];
	const mxArray * output_filename_mx = prhs[2];
  
	const int * size_F = mxGetDimensions(F_mx);
	int n_node = size_F[0];
	double * F = mxGetPr(F_mx);

  double * adjacency_0 = mxGetPr(adjacency_mx);
  double * adjacency = (double *) new double[n_node*n_node];
  for(int i=0; i<n_node; i++)
    for(int j=0; j<n_node; j++)
      adjacency[i+ j*n_node] = adjacency_0[i + j*n_node];
  /*
  for(int i=0; i<n_node; i++)
    for(int j=0; j<n_node; j++){
      if(adjacency_0[i + j*n_node]==0)
        continue;
      for(int k=0; k<n_node; k++)
        if(adjacency_0[j + k*n_node]!=0)
          adjacency[i + k*n_node] = adjacency[k + i*n_node] = 1;
    }
	*/
	char * output_filename = mxArrayToString(output_filename_mx);
  mexPrintf("output_filename: %s\n", output_filename);
  
  FILE * fout = fopen(output_filename, "wt");
  fprintf(fout, "min: ");
  
  for(int i=0; i<n_node; i++){
    for(int j=i+1; j<n_node; j++){
      if(adjacency[i + j*n_node]==0)
        continue;
      fprintf(fout, "+ %g x_%d_%d ", F[i + j*n_node], i, j);
    }
  }
  fprintf(fout, ";\n");

  for(int i=0; i<n_node; i++){
    for(int j=i+1; j<n_node; j++){
      if(adjacency[i + j*n_node]==0)
        continue;
      fprintf(fout, "0 <= x_%d_%d <= 1;\n", i, j);
    }
  }

  for(int i=0; i<n_node; i++){
    for(int j=i+1; j<n_node; j++){
      if(adjacency[i + j*n_node]==0)
        continue;
      for(int k=0; k<i; k++){
        if(adjacency[k + i*n_node]==0 || adjacency[k + j*n_node]==0)
          continue;
        fprintf(fout, "x_%d_%d <= x_%d_%d + x_%d_%d;\n", i, j, k, i, k, j);
        fprintf(fout, "x_%d_%d <= x_%d_%d + x_%d_%d;\n", k, i, i, j, k, j);
        fprintf(fout, "x_%d_%d <= x_%d_%d + x_%d_%d;\n", k, j, k, i, i, j);
      }
      for(int k=i+1; k<j; k++){
        if(adjacency[i + k*n_node]==0 || adjacency[k + j*n_node]==0)
          continue;
        fprintf(fout, "x_%d_%d <= x_%d_%d + x_%d_%d;\n", i, j, i, k, k, j);
        fprintf(fout, "x_%d_%d <= x_%d_%d + x_%d_%d;\n", i, k, i, j, k, j);
        fprintf(fout, "x_%d_%d <= x_%d_%d + x_%d_%d;\n", k, j, i, k, i, j);
      }
      for(int k=j+1; k<n_node; k++){
        if(adjacency[i + k*n_node]==0 || adjacency[j + k*n_node]==0)
          continue;
        fprintf(fout, "x_%d_%d <= x_%d_%d + x_%d_%d;\n", i, j, i, k, j, k);
        fprintf(fout, "x_%d_%d <= x_%d_%d + x_%d_%d;\n", i, k, i, j, j, k);
        fprintf(fout, "x_%d_%d <= x_%d_%d + x_%d_%d;\n", j, k, i, k, i, j);
      }
    }
  }

  fclose(fout);
  
  delete [] adjacency;
	mxFree(output_filename);
  
  return;
}
