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
      plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
      *mxGetPr(plhs[0]) = 1;
      return;
    }
    mexPrintf("Usage segment_map = remove_merged_boundaries_2D(segment_map)\n");
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
  
  //mexPrintf("START: remove_merged_boundaries_2D\n");
  int numDim = mxGetNumberOfDimensions(prhs[0]);
  if(numDim!=2){
    mexErrMsgTxt("Wrong no. of dimensions for arg. 1\n");
    return;
  }
  
  const int * sizeImage;
  sizeImage = mxGetDimensions(prhs[0]);
  int width = sizeImage[1], height = sizeImage[0];
  int n_pixel = height*width;
  unsigned int * segment_map = (unsigned int *) mxGetPr(prhs[0]);
  plhs[0] = mxCreateNumericMatrix(height, width, mxUINT32_CLASS, mxREAL);
  unsigned int * output_label_map = (unsigned int *) mxGetPr(plhs[0]);
  
  {
    int x, y, xmin, xmax, ymin, ymax, x1, y1;
    int z_o, x_o, y_o;
    unsigned int *l, *o, s, n[4], nn, ni;
    bool is_border, is_set;
    l = segment_map;
    o = output_label_map;
    for(x=0; x<width; x++)
      for(y=0; y<height; y++, l++, o++)
        *o = *l;
    
    l = segment_map;
    o = output_label_map;
    for(x=0; x<width; x++){
      for(y=0; y<height; y++, l++, o++){
        if(*l!=0){
          continue;
        }

        is_border = false;
        if(x>0 && x<width-1 && *(l-height)!=0 && *(l+height)!=0 
                && *(l-height)!=*(l+height))
          is_border = true;
        else{
          if(y>0 && y<height-1 && *(l-1)!=0 && *(l+1)!=0
                  && *(l-1)!=*(l+1))
            is_border = true;
          else{
            if(x>0 && x<width-1 && y>0 && y<height-1 && 
                    *(l-height)!=0 && *(l+height)!=0 && *(l-1)!=0 && *(l+1)!=0 &&
                    *(l-height)==*(l+height) && *(l-1)==*(l+1) &&
                    *(l-height)!=*(l-1))
              is_border = true;
          }
        }
        
        if(is_border)
          *o = 0;
        else{
          if(x>0 && x<width-1 && *(l-height)!=0 && *(l+height)!=0
                  && *(l-height)==*(l+height))
            *o = *(l-height);
          if(y>0 && y<height-1 && *(l-1)!=0 && *(l+1)!=0
                  && *(l-1)==*(l+1))
            *o = *(l-1);
        }
      }
    }
  }
  
  //mexPrintf("STOP: remove_merged_boundaries_2D\n");
  return;
}
