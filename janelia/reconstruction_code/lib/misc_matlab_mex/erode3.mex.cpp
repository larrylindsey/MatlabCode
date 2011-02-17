// compute 3D median filter
// Shiv N. Vitaladevuni
//

#include <mex.h>
#include <math.h>
#include <algorithm>

#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

/*
 * input params:
 *	1.	PxQxR 3D matrix
 *	2.	v - size of median filter kernel (2v+1)x(2v+1)x(2v+1)
 *
 * output:
 *	Computes pxq 2D histogram of the vectors
 *
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    mexPrintf("input params:\n");
		mexPrintf("\t1.	PxQxR 3D matrix\n");
		mexPrintf("\t2.	v - size of median filter kernel (2v+1)x(2v+1)x(2v+1)\n");
    mexPrintf("output:\n");
		mexPrintf("\tComputes 3D median filter as a PxQxR 3D matrix\n");
    return;
  }
  if(nrhs!=2)
  {
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

	const mxArray * volume_mx = prhs[0];
	const mxArray * kernel_size_mx = prhs[1];
  
	const int * size_volume = mxGetDimensions(volume_mx);
	int height = size_volume[0];
	int width = size_volume[1];
	int depth = size_volume[2];
	double * volume = mxGetPr(volume_mx);
	
	int size_kernel = (int) (*mxGetPr(kernel_size_mx));

  plhs[0] = mxCreateNumericArray(3, size_volume, mxDOUBLE_CLASS, mxREAL);
  double * median_filter_p = mxGetPr(plhs[0]);
  
  double min_value;
  
  int n_pixel_per_plane = height*width;
  int y_skip;
  int plane_skip;
  double * v_p;
  for(int d=0; d<depth; d++){
    for(int x=0; x<width; x++){
      for(int y=0; y<height; y++, median_filter_p++){
        int d0 = MAX(0, d-size_kernel);
        int d1 = MIN(depth-1, d+size_kernel);
        int x0 = MAX(0, x-size_kernel);
        int x1 = MIN(width-1, x+size_kernel);
        int y0 = MAX(0, y-size_kernel);
        int y1 = MIN(height-1, y+size_kernel);
        
        y_skip = height - y1 + y0 - 1;
        plane_skip = n_pixel_per_plane - (x1-x0)*height - y1 + y0 - 1 - y_skip;
        
        v_p = volume + d0*n_pixel_per_plane + x0*height + y0;
        
        min_value = *v_p;
        
        for(int dk=d0; dk<=d1; dk++, v_p+=plane_skip){
          for(int xk=x0; xk<=x1; xk++, v_p+=y_skip){
            for(int yk=y0; yk<=y1; yk++, v_p++){
              min_value = MIN(min_value, *v_p);
            }
          }
        }
        
        *median_filter_p = min_value;
      }
    }
  }
  
  return;
}
