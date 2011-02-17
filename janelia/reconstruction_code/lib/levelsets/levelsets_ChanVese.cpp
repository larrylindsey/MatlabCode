#include "stdio.h"
#include "math.h"
#include "stdlib.h"

#include <levelsets_ChanVese.h>

# define eps 2.2204e-16
# define PI 3.1416

void ChanVese(double *phi, double *phiNew, double *Im, double epsilon, double nu, double delt, int r, int c, double c1, double c2, double lambda1, double lambda2, int h, int n_iteration) {
  
  double DataTerm, CurvatureTerm;
  double **phiShiftedLeft, **phiShiftedRight, **phiShiftedUp, **phiShiftedDown, deleps;
  double C, C1, C2, C3, C4, m;
  double **phiForwardx, **phiForwardy, **phiBackwardx, **phiBackwardy;
  double **phiCentralx, **phiCentraly, **phiCentralyShiftedRight, **phiCentralxShiftedDown;
  int i, j;
  
  for (i=0; i<r; i++) {
    for (j=0; j<c; j++) {
      phiNew[j*r+i] = phi[j*r+i];
    }
  }
  
  phiForwardx=new double *[r];
  phiForwardy=new double *[r];
  phiBackwardx=new double *[r];
  phiBackwardy=new double *[r];
  phiCentralx=new double *[r];
  phiCentraly=new double *[r];
  phiCentralyShiftedRight=new double *[r];
  phiCentralxShiftedDown=new double *[r];
  phiShiftedLeft=new double *[r];
  phiShiftedRight=new double *[r];
  phiShiftedUp=new double *[r];
  phiShiftedDown=new double *[r];
  
  for (i=0; i<r; i++) {
    phiForwardx[i]=new double [c];
    phiForwardy[i]=new double [c];
    phiBackwardx[i]=new double [c];
    phiBackwardy[i]=new double [c];
    phiCentralx[i]=new double[c];
    phiCentraly[i]=new double[c];
    phiCentralyShiftedRight[i]=new double[c];
    phiCentralxShiftedDown[i]=new double [c];
    phiShiftedLeft[i]=new double [c];
    phiShiftedRight[i]=new double [c];
    phiShiftedUp[i]=new double [c];
    phiShiftedDown[i]=new double [c];
  }

  int iter;
  for(iter=0; iter<n_iteration; iter++){
    for (i=0; i<r; i++) {
      for (j=0; j<c; j++) {
        if (j==(c-1)) {
          phiForwardx[i][j]=0.0;
          phiShiftedLeft[i][j]=phiNew[(j*r)+i];
        }
        else {
          phiForwardx[i][j]=phiNew[((j+1)*r)+i]-phiNew[(j*r)+i];
          phiShiftedLeft[i][j]=phiNew[((j+1)*r)+i];
        }
      }
    }
  
    for (i=0; i<r; i++) {
      for (j=0; j<c; j++) {
        if (j==0) {
          phiBackwardx[i][j]=0.0;
          phiShiftedRight[i][j]=phiNew[(j*r)+i];
        }
        else {
          phiBackwardx[i][j]=phiNew[(j*r)+i]-phiNew[((j-1)*r)+i];
          phiShiftedRight[i][j]=phiNew[((j-1)*r)+i];
        }
      }
    }
  
  
    for (i=0; i<r; i++) {
      for (j=0; j<c; j++) {
        if (i==(r-1)) {
          phiForwardy[i][j]=0.0;
          phiShiftedUp[i][j]=phiNew[(j*r)+i];
        }
        else {
          phiForwardy[i][j]=phiNew[(j*r)+i+1]-phiNew[(j*r)+i];
          phiShiftedUp[i][j]=phiNew[(j*r)+i+1];
        }
      }
    }
  
    for (i=0; i<r; i++) {
      for (j=0; j<c; j++) {
        if (i==0) {
          phiBackwardy[i][j]=0.0;
          phiShiftedDown[i][j]=phiNew[(j*r)+i];
        }
        else {
          phiBackwardy[i][j]=phiNew[(j*r)+i]-phiNew[(j*r)+i-1];
          phiShiftedDown[i][j]=phiNew[(j*r)+i-1];
        }
      }
    }
  
    for (i=0; i<r; i++) {
      for (j=0; j<c; j++) {
        if (j==0) {
          phiCentralx[i][j]=phiNew[((j+1)*r)+i]-phiNew[(j*r)+i];
        }
        else if (j==(c-1)) {
          phiCentralx[i][j]=phiNew[(j*r)+i]-phiNew[((j-1)*r)+i];
        }
        else {
          phiCentralx[i][j]=phiNew[((j+1)*r)+i]-phiNew[((j-1)*r)+i];
        }
      }
    }
  
    for (i=0; i<r; i++) {
      for (j=0; j<c; j++) {
        if (i==0) {
          phiCentraly[i][j]=phiNew[(j*r)+i+1]-phiNew[(j*r)+i];
        }
        else if (i==(r-1)) {
          phiCentraly[i][j]=phiNew[(j*r)+i]-phiNew[(j*r)+i-1];
        }
        else {
          phiCentraly[i][j]=phiNew[(j*r)+i+1]-phiNew[(j*r)+i-1];
        }
      }
    }
  
    for (i=0; i<r; i++) {
      for (j=0; j<c; j++) {
        if(j==0)
          phiCentralyShiftedRight[i][j]=phiCentraly[i][j];
        else
          phiCentralyShiftedRight[i][j]=phiCentraly[i][j-1];
      }
    }
  
    for (i=0; i<r; i++) {
      for (j=0; j<c; j++) {
        if(i==0)
          phiCentralxShiftedDown[i][j]=phiCentralx[i][j];
        else
          phiCentralxShiftedDown[i][j]=phiCentralx[i-1][j];
      }
    }
  
    for (i=0; i<r; i++) {
      for (j=0; j<c; j++) {
        C1=1/(sqrt(pow(phiForwardx[i][j]/h, 2)+pow(phiCentraly[i][j]/(2*h), 2))+eps);
        C2=1/(sqrt(pow(phiBackwardx[i][j]/h, 2)+pow(phiCentralyShiftedRight[i][j]/(2*h), 2))+eps);
        C3=1/(sqrt(pow(phiCentralx[i][j]/(2*h), 2)+pow(phiForwardy[i][j]/h, 2))+eps);
        C4=1/(sqrt(pow(phiCentralxShiftedDown[i][j]/(2*h), 2)+pow(phiBackwardy[i][j]/h, 2))+eps);
        
        deleps=epsilon/(PI*(pow(epsilon, 2)+pow(phiNew[(j*r)+i], 2)));
        m=delt*nu*deleps/(h*h);
        C=1.0+m*(C1+C2+C3+C4);
      
        DataTerm=-lambda1*(Im[(j*r)+i]-c1)*(Im[(j*r)+i]-c1)+lambda2*(Im[(j*r)+i]-c2)*(Im[(j*r)+i]-c2);
        CurvatureTerm=m*(C1*phiShiftedLeft[i][j]+C2*phiShiftedRight[i][j]+C3*phiShiftedUp[i][j]+C4*phiShiftedDown[i][j]);
        phiNew[(j*r)+i]=(phiNew[(j*r)+i]+CurvatureTerm+delt*deleps*DataTerm)/C;
      }
    }
  }

  
  
  for (i=0; i<r; i++) {
    delete []  phiForwardx[i];
    delete []  phiForwardy[i];
    delete []  phiBackwardx[i];
    delete []  phiBackwardy[i];
    delete []  phiCentralx[i];
    delete []  phiCentraly[i];
    delete []  phiCentralyShiftedRight[i];
    delete []  phiCentralxShiftedDown[i];
    delete []  phiShiftedLeft[i];
    delete []  phiShiftedRight[i];
    delete []  phiShiftedUp[i];
    delete []  phiShiftedDown[i];
  }
  
  delete [] phiForwardx;
  delete [] phiForwardy;
  delete [] phiBackwardx;
  delete [] phiBackwardy;
  delete [] phiCentralx;
  delete [] phiCentraly;
  delete [] phiCentralyShiftedRight;
  delete [] phiCentralxShiftedDown;
  delete [] phiShiftedLeft;
  delete [] phiShiftedRight;
  delete [] phiShiftedUp;
  delete [] phiShiftedDown;
  
}
        
