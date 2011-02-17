#include "stdio.h"
#include "math.h"
# include "stdlib.h"

#include <levelsets_ChanVese.h>

void foo(void) {
  
  double *phi, *Im;
  double epsilon=1.0, nu, delt=0.1, c1=0.0, c2=255.0, lambda1, lambda2;
  double *phiNew, *phiTemp;
  int i, j, k, r, c, h=1, ITER;
  
  
  ChanVese(phiTemp, phiNew, Im, epsilon, nu, delt, r, c, c1, c2, lambda1, lambda2, h);
    

  return;
}
        
