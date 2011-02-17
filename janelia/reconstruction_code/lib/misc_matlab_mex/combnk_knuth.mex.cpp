// generate k combinations of n numbers.
// Knuth's algorithm
//

#include <mex.h>
#include <stdlib.h> 


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
		mexPrintf("\t1.	n\n");
		mexPrintf("\t2.	k\n");
    mexPrintf("output:\n");
		mexPrintf("\tCombinations of n\n");
    return;
  }
  if(nrhs!=2)
  {
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

	int i, j=1, k, n, *c, x;
  n = (int) * mxGetPr(prhs[0]);
  k = (int) * mxGetPr(prhs[1]);
	c = (int *) new int[k+3];
  
  int n_comb=1;
  for(i=n; i>n-k; i--)
    n_comb *= i;
  for(i=k;i>0; i--)
    n_comb /= i;
  plhs[0] = mxCreateNumericMatrix(n_comb, k, mxUINT8_CLASS, mxREAL);
  unsigned char * comb = (unsigned char *) mxGetPr(plhs[0]);
  int comb_i=0;

  
  for (i=1; i <= k; i++) c[i] = i;
	c[k+1] = n+1;
	c[k+2] = 0;
	j = k;
  visit:
	for (i=k; i >= 1; i--)
    comb[comb_i+(i-1)*n_comb] = c[i];
	comb_i++;
	if (j > 0) {x = j+1; goto incr;}
	if (c[1] + 1 < c[2])
  {
    c[1] += 1;
    goto visit;
  }
	j = 2;
  do_more:
	c[j-1] = j-1;
	x = c[j] + 1;
	if (x == c[j+1]) {j++; goto do_more;}
	if (j > k) return;
  incr:
	c[j] = x;
	j--;
	goto visit;

  delete [] c;
  
	return;
}
