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

#define MX_CLASS mxUINT8_CLASS
typedef  unsigned char Data_Type;

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
		mexPrintf("\t2.	[u,v,w] - size of median filter (2u+1)x(2v+1)x(2w+1)\n");
    mexPrintf("output:\n");
		mexPrintf("\tComputes 3D median filter as a PxQxR 3D matrix\n");
    return;
  }
  if(nrhs!=2)
  {
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

  mexPrintf("START: medfilt3\n");
	const mxArray * volume_mx = prhs[0];
	const mxArray * kernel_size_mx = prhs[1];
  
	const int * size_volume = mxGetDimensions(volume_mx);
	int height = size_volume[0];
	int width = size_volume[1];
	int depth = size_volume[2];
	Data_Type * volume = (Data_Type *) mxGetPr(volume_mx);
	
	double * size_kernel = mxGetPr(kernel_size_mx);
  int size_kernel_u = size_kernel[0];
  int size_kernel_v = size_kernel[1];
  int size_kernel_w = size_kernel[2];
  mexPrintf("size_kernel_u: %d, size_kernel_v: %d, size_kernel_w: %d\n",
            size_kernel_u, size_kernel_v, size_kernel_w);
  
  plhs[0] = mxCreateNumericArray(3, size_volume, MX_CLASS, mxREAL);
  Data_Type * median_filter_p = (Data_Type *) mxGetPr(plhs[0]);
  
  Data_Type * neigh_values = (Data_Type *) new Data_Type[(2*size_kernel_u+1)*(2*size_kernel_v+1)*(2*size_kernel_w+1)];
  
  int n_pixel_per_plane = height*width;
  
  for(int d=0; d<depth; d++){
    mexPrintf("d: %d\n", d);
    mexEvalString("drawnow;");
    for(int x=0; x<width; x++){
      for(int y=0; y<height; y++, median_filter_p++){
        int d0 = MAX(0, d-size_kernel_u);
        int d1 = MIN(depth-1, d+size_kernel_u);
        int x0 = MAX(0, x-size_kernel_v);
        int x1 = MIN(width-1, x+size_kernel_v);
        int y0 = MAX(0, y-size_kernel_w);
        int y1 = MIN(height-1, y+size_kernel_w);
        
        int y_skip = height - y1 + y0 - 1;
        int plane_skip = n_pixel_per_plane - (x1-x0)*height -
          y1 + y0 - 1 - y_skip;
        
        Data_Type * v_p = volume + d0*n_pixel_per_plane + x0*height + y0;
        Data_Type * n_v_p = neigh_values;
        int n_neigh = 0;
        for(int dk=d0; dk<=d1; dk++, v_p+=plane_skip){
          for(int xk=x0; xk<=x1; xk++, v_p+=y_skip){
            for(int yk=y0; yk<=y1; yk++, v_p++, *n_v_p++, n_neigh++){
              *n_v_p = *v_p;
            }
          }
        }
        
        std::nth_element(neigh_values, neigh_values+n_neigh/2,
                         neigh_values+n_neigh-1);
        *median_filter_p = neigh_values[n_neigh/2];
      }
    }
  }

  delete [] neigh_values;
  mexPrintf("STOP: medfilt3\n");
  return;
}
