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
 *	1.	image intensity map height x width
 *	2.	sampling mask height x width
 *	3.	neighborhood window sizes - N x 1 vector
 *	4.	intensity bins - M x 1 vector
 *
 * output:
 *	Computes M*N dimensional feature vectors of the intensity histograms
 *	1. R X M*N matrix, where R is the number of samples in the input sampling mask.
 *	2. 2 X M*N matrix, coordinates of the samples on the image plane.
 *
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    mexPrintf("input params:\n");
    mexPrintf("\t1.	image intensity map height x width\n");
    mexPrintf("\t2.	sampling mask height x width\n");
    mexPrintf("\t3.	neighborhood window sizes - N x 1 vector\n");
    mexPrintf("\t4.	intensity bins - M x 1 vector\n\n");
    mexPrintf("output:\n");
    mexPrintf("\tComputes M*N dimensional feature vectors of the intensity histograms\n");
    mexPrintf("\t1. R X M*N matrix, where R is the number of samples in the input sampling mask.\n");
    mexPrintf("\t2. R X 2 matrix, coordinates of the samples on the image plane.\n");
    return;
  }
  if(nrhs!=4)
  {
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

  if(nlhs!=2)
  {
    mexErrMsgTxt("Wrong number of outputs\n");
    return;
  }

	const mxArray * intensity_mx = prhs[0];
  const mxArray * sampling_mask_mx = prhs[1];
  const mxArray * window_sizes_mx = prhs[2];
	const mxArray * intensity_bins_mx = prhs[3];
  
  const int * size_image = mxGetDimensions(intensity_mx);
  int height = size_image[0];
	int width = size_image[1];
	int n_pixel = height*width;
	double * intensity_map = mxGetPr(intensity_mx);
	
	const int * size_neigh_windows = mxGetDimensions(window_sizes_mx);
	int n_window = size_neigh_windows[0]*size_neigh_windows[1];
	double * window_sizes = mxGetPr(window_sizes_mx);
	int max_window_size=-1;
  double * window_areas = (double *) new double[n_window];
	for(int i=0; i<n_window; i++){
		max_window_size = (int) MAX(max_window_size,window_sizes[i]);
    window_areas[i] = (2*window_sizes[i]+1)*(2*window_sizes[i]+1);
	}
	  
	const int * size_intensity_bins = mxGetDimensions(intensity_bins_mx);
	int n_intensity_bin = size_intensity_bins[0]*size_intensity_bins[1];
	double * intensity_bins = mxGetPr(intensity_bins_mx);
	double * bin_boundaries = (double *) new double[n_intensity_bin+1];
	bin_boundaries[0]=0.0;
	for(int i=1; i<n_intensity_bin; i++){
		bin_boundaries[i] = (intensity_bins[i-1]+intensity_bins[i])/2.0;
	}
	bin_boundaries[n_intensity_bin] = 1.0;
	
	double * sample = mxGetPr(sampling_mask_mx);
	double * sample_mask = (double *) new double[n_pixel];
	for(int i=0; i<n_pixel; i++)
		sample_mask[i] = sample[i];
	
	//remove samples near the border
	for(int y=0; y<height; y++)
		for(int x=0; x<max_window_size+5; x++)
			sample_mask[(width-x-1)*height+y]=sample_mask[x*height+y]=0;

	
	for(int x=0; x<width; x++)
		for(int y=0; y<max_window_size+5; y++)
			sample_mask[x*height+y] = sample_mask[x*height+height-y-1]=0;
	
	
	int n_sample = 0;
	double * s_p = sample_mask;
	for(int x=0; x<width; x++)
		for(int y=0; y<height; y++, s_p++)
			if(*s_p==1)
				n_sample++;
	plhs[0] = mxCreateDoubleMatrix(n_sample, n_window*n_intensity_bin, mxREAL);
	double * feature = mxGetPr(plhs[0]);
	
	plhs[1] = mxCreateDoubleMatrix(n_sample, 2, mxREAL);
	double * sample_loc = mxGetPr(plhs[1]);
	s_p = sample_mask;
	double * s_l_p = sample_loc;
	for(int x=0; x<width; x++)
		for(int y=0; y<height; y++, s_p++)
			if(*s_p==1){
				*s_l_p = x;
				*(s_l_p+n_sample)=y;
				s_l_p++;
			}
	//
	//compute integral images
	//
	double * cum_hist = (double *) new double[n_pixel];
	double * i_p, * c_p;
	double curr_sum;
	double * f_p = feature;
	double * c_p_1, * c_p_2, * c_p_3, * c_p_4;
	for(int bin=0; bin<n_intensity_bin; bin++){
		cum_hist[0]=(intensity_map[0]>=bin_boundaries[bin] && intensity_map[0]<bin_boundaries[bin+1])?1:0;
			
		i_p = intensity_map+1;
		c_p = cum_hist+1;
		for(int y=1; y<height; y++, i_p++, c_p++){
			*c_p = *(c_p-1) + ((*i_p>=bin_boundaries[bin] && *i_p<bin_boundaries[bin+1])?1:0);
		}
		for(int x=1; x<width; x++){
			cum_hist[x*height]=cum_hist[(x-1)*height] + ((intensity_map[x*height]>=bin_boundaries[bin] && intensity_map[x*height]<bin_boundaries[bin+1])?1:0);
				
			curr_sum = cum_hist[x*height];
				
			i_p++;
			c_p++;
			for(int y=1; y<height; y++, i_p++, c_p++){
				curr_sum += (*i_p>=bin_boundaries[bin] && *i_p<bin_boundaries[bin+1])?1:0;
				*c_p = *(c_p-height) + curr_sum;
			}
		}
		
		for(int w=0; w<n_window; w++){
			
			s_p = sample_mask;
			c_p_1 = cum_hist + (int)(-(window_sizes[w]+1)*height - window_sizes[w]-1);
			c_p_2 = cum_hist + (int)(window_sizes[w]*height + window_sizes[w]);
			c_p_3 = cum_hist + (int)(-(window_sizes[w]+1)*height + window_sizes[w]);
			c_p_4 = cum_hist + (int)(window_sizes[w]*height - window_sizes[w] -1);
			for(int x=0; x<width; x++){
				for(int y=0; y<height; y++, s_p++, c_p_1++, c_p_2++, c_p_3++, c_p_4++){
					if(*s_p==1){
						*f_p = (*c_p_1 + *c_p_2 - *c_p_3 - *c_p_4)/window_areas[w];
						f_p++;
					}
				}
			}
		}
	}

  delete [] window_areas;
	delete [] cum_hist;
	delete [] sample_mask;
  delete [] bin_boundaries;
  
  return;
}
