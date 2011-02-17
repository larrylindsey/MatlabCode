#include "stdio.h"
#include "math.h"
# include "stdlib.h"

#include "mex.h"  

#include <levelsets_ChanVese.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  if(nrhs==0){
    mexPrintf("Usage: phi_new = levelsets_ChanVese(blah blah blah);\n");
    return;
  }
  double *phi, *Im;
  double epsilon=1.0, nu, delt=0.1, c1=0.0, c2=255.0, lambda1, lambda2;
  double *phiNew, *phiTemp;
  int i, j, k, r, c, h=1, ITER;
  
  
  nu=mxGetScalar(prhs[0]);
  phi = mxGetPr(prhs[1]);
  r = (int)mxGetM(prhs[1]);
  c = (int)mxGetN(prhs[1]);
  Im = mxGetPr(prhs[2]);
  lambda1=mxGetScalar(prhs[3]);
  lambda2=mxGetScalar(prhs[4]);
  ITER=(int)mxGetScalar(prhs[5]);
  
  
  plhs[0] = mxCreateDoubleMatrix(r, c, mxREAL);
  phiNew=mxGetPr(plhs[0]);
  phiTemp=new double [r*c];
  
 
  
  for (i=0; i<r; i++)
    for (j=0; j<c; j++)
      phiTemp[(j*r)+i]=-phi[(j*r)+i];
  
  for (i=0;i<ITER;i++) {
    ChanVese(phiTemp, phiNew, Im, epsilon, nu, delt, r, c, c1, c2, lambda1, lambda2, h);
    
    for (j=0; j<r; j++)
      for (k=0; k<c; k++)
        phiTemp[(k*r)+j]=phiNew[(k*r)+j];
    
    
  }
  
  for (i=0; i<r; i++)
    for (j=0; j<c; j++)
      phiNew[(j*r)+i]=-phiNew[(j*r)+i]; 
  
  
  delete [] phiTemp;
  return;
}
        
