// Remove boundaries (0 label pixels) between merged segments
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	03312009	init. code
//

#include <mex.h>
#include <vector>

#define __MATLAB_MEX__
#include <image_2D_utilities.h>


void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    if(nlhs==1){
      plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
      *mxGetPr(plhs[0]) = 1;
      return;
    }
    mexPrintf("Usage segment_map = remove_merged_boundaries_3D(segment_map)\n");
    mexPrintf("input params:\n");
    mexPrintf("\t1. segment label map MxN matrix (uint32)\n");
    mexPrintf("output:\n");
    mexPrintf("\t1. segment label map with merged boundaries erased MxN matrix (uint32)\n");
    return;
  }
  if(nrhs!=1){
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }
  if(nlhs>1){
    mexErrMsgTxt("Wrong number of outputs\n");
    return;
  }

  //mexPrintf("START: remove_merged_boundaries_3D\n");
  int numDim = mxGetNumberOfDimensions(prhs[0]);
  if(numDim!=3){
    mexErrMsgTxt("Wrong no. of dimensions for arg. 1\n");
    return;
  }

  const int * sizeImage;
  sizeImage = mxGetDimensions(prhs[0]);
  int width = sizeImage[1], height = sizeImage[0], depth = sizeImage[2];
  int n_pixel = height*width*depth;
  int n_pixel_plane = height*width;
  unsigned int * segment_map = (unsigned int *) mxGetPr(prhs[0]);
  plhs[0] = mxCreateNumericArray(3, sizeImage, mxUINT32_CLASS, mxREAL);
  unsigned int * output_label_map = (unsigned int *) mxGetPr(plhs[0]);

  {
    int x, y, z, xmin, xmax, ymin, ymax, zmin, zmax, x1, y1, z1;
    int z_o, x_o, y_o;
    unsigned int *l, *o, s;
    bool is_border;
    l = segment_map;
    o = output_label_map;
    for(z=0; z<depth; z++){
      for(x=0; x<width; x++){
	for(y=0; y<height; y++, l++, o++){
	  if(*l!=0){
	    *o=*l;
	    continue;
	  }

	  xmin = std::max(0, x-1);
	  xmax = std::min(width-1, x+1);
	  ymin = std::max(0, y-1);
	  ymax = std::min(height-1, y+1);
	  zmin = std::max(0, z-1);
	  zmax = std::min(depth-1, z+1);

	  s = 0;
	  is_border = false;
	  for(z1=zmin; z1<=zmax; z1++){
	    z_o = z1*n_pixel_plane;
	    for(x1=xmin; x1<=xmax; x1++){
	      x_o = z_o + x1*height;
	      for(y1=ymin; y1<=ymax; y1++){
		y_o = x_o + y1;
		if(segment_map[y_o]!=0){
		  if(s!=0 && s!=segment_map[y_o]){
		    is_border = true;
		    break;
		  }
		  else{
		    if(s==0)
		      s = segment_map[y_o];
		  }
		}
	      }
	      if(is_border)
		break;
	    }
	    if(is_border)
	      break;
	  }
	  if(is_border)
	    *o = 0;
	  else
	    *o = s;
	}
      }
    }
  }
  
  //mexPrintf("STOP: remove_merged_boundaries_3D\n");
  return;
}
