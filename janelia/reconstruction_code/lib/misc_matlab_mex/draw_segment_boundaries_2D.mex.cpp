// draw segment boundaries for a given segmentation label map
// assumes that the segments are tight, i.e. no gaps

#include <mex.h>
#include <math.h>

#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

#define MIN(A,B)	(A<B?A:B)
#define MAX(A,B)	(A>B?A:B)

#define NUM_ORIENT_BINS		180

/*
 * input params:
 *	label map - 2D matrix
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0 && nlhs==0){
    mexPrintf("label_map_with_boundaries = draw_segment_boundaries_2D(label_map)\n");
    mexPrintf("Input:\n");
    mexPrintf("\t1: segmentation label map with unique labels (uint32 MxN matrix).\n");
    mexPrintf("Output:\n");
    mexPrintf("\t1: segmentation label map (uint32 MxN matrix) with boundaries (label=0)\n");
    mexPrintf("\t\tdrawn along the segment edges.\n");
    return;
  }
  
	if(nrhs!=1){
		mexErrMsgTxt("Wrong number of inputs\n");
		return;
	}


	int numDim = mxGetNumberOfDimensions(prhs[0]);
	if(numDim!=2){
		mexErrMsgTxt("Wrong no. of dimensions for arg. 1\n");
		return;
	}

	const int * sizeImage;
	sizeImage = mxGetDimensions(prhs[0]);
	int width = sizeImage[1], height = sizeImage[0];

	unsigned int * label_map = (unsigned int *) mxGetPr(prhs[0]);

	int x, y;

	plhs[0] = mxCreateNumericMatrix(height, width, mxUINT32_CLASS, mxREAL);
	unsigned int * boundary_map = (unsigned int *) mxGetPr(plhs[0]);
	unsigned int * boundary_map_ptr = boundary_map;
	unsigned int  * label_map_ptr = label_map;
	for(x=0; x<width-1; x++, label_map_ptr++,boundary_map_ptr++){
		for(y=0; y<height-1; y++, label_map_ptr++,boundary_map_ptr++){
      if( (*label_map_ptr)==0){
        *boundary_map_ptr=0;
        continue;
      }
			if( (*label_map_ptr)>0 && *(label_map_ptr+1)>0 && ((*label_map_ptr)!=(*(label_map_ptr+1)))){
        *boundary_map_ptr=0;
        continue;
      }
			if( (*label_map_ptr)>0 && *(label_map_ptr+height)>0 && ((*label_map_ptr)!=(*(label_map_ptr+height)))){
        *boundary_map_ptr=0;
        continue;
      }
  		*boundary_map_ptr=*label_map_ptr;
		}
	}

	for(x=0, boundary_map_ptr=boundary_map+height-1, label_map_ptr=label_map+height-1; x<width-1; x++, boundary_map_ptr+=height, label_map_ptr+=height){
		if( (*label_map_ptr)>0 && (*(label_map_ptr+height))>0 && (*label_map_ptr)!=(*(label_map_ptr+height)) )
			*boundary_map_ptr=0;
		else
			*boundary_map_ptr=*label_map_ptr;
	}

	for(y=0, boundary_map_ptr=boundary_map+(width-1)*height, label_map_ptr=label_map+(width-1)*height; y<height-1; y++, boundary_map_ptr++, label_map_ptr++){
		if( (*label_map_ptr)>0 && (*(label_map_ptr+1))>0 && (*label_map_ptr)!=(*(label_map_ptr+1)) )
			*boundary_map_ptr=0;
		else
			*boundary_map_ptr=*label_map_ptr;
	}

  boundary_map[width*height-1] = label_map[width*height-1];
	return;

}
