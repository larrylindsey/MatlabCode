/*
 * Code for 3D watershed constrained by a sub-volume mask and seeded in user-defined
 * z-plane. Called from proofreading gui.
 *
 * Uses a modification of Gene Myers' (JFRC) 3D watershed code. For license contact
 * Dr. Gene Myers, JFRC, HHMI.
 *
 * v0   11122008  init code for testing
 * v1   11132008  code modified for proofreading gui
 *
 * Shiv N. Vitaladevuni
 * Janelia Farm Research Campus, HHMI.
 *
 */

#include <mex.h>
#include <math.h>
//#define inline 
//#define const 
#include <image_lib.h>
#include <water_shed.h>

#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

#define MIN(A,B)	(A<B?A:B)
#define MAX(A,B)	(A>B?A:B)

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  const mxArray * mode_mx, * image_stack_mx, * image_stack_scale_mx, * cat_mx, * superpixel_2_seg_map_mx, * pmap_mx, * mask_label_mx, * seed_image_mx, * seed_depth_mx, * label_new_stack_mx;
  int mode;
  const int * size_0;
  int height, width, depth;
  Stack * image_stack;
  int image_stack_scale;
  int i;
  unsigned short int ** superpixel_maps;
  const int * size_1;
  int original_width;
  double ** superpixel_2_seg_maps;
  unsigned int * pmap;
  unsigned int mask_label;
  Image * seed_image;
  unsigned int seed_depth;
  Stack * label_new;
  
	if(nrhs==0){
		mexPrintf("Usage gui_3D_split_watershed(mode, image_stack, cat, superpixel_2_seg_map, pmap, ...\n");
    mexPrintf("\tpmap, mask_label, seed_image, seed_depth. label_new_stack);\n");
		mexPrintf("input params:\n");
		mexPrintf("\t1. mode: 0 initialize, 1 segment, 2 clear data structures.\n");
		mexPrintf("\t2. image_stack: 3D stack of intensity volume MxNxP uint8 matrix\n");
		mexPrintf("\t3. image_stack_scale (s): down-scaling along x and y axis\n");
		mexPrintf("\t4. cat{}: cell array of P 2D label maps volume sMxsN uint16 matrices\n");
		mexPrintf("\t5. superpixel_2_seg_map: cell array of superpixel to segment mapping\n");
		mexPrintf("\t6. pmap: segment to body id\n");
		mexPrintf("\t7. mask_label: body id to within which segmentation is to be performed uint32\n");
		mexPrintf("\t8. seed_image: image of connected components as seeds MxN uint8\n");
		mexPrintf("\t9. seed_depth: depth of seed image in the stack uint32\n");
		mexPrintf("\t10. label_new_stack: 3D stack of new label volume MxNxP int32 matrix\n");
		mexPrintf("output:\n");
		mexPrintf("\t1. label_stack_new: 3D stack of label volume MxNxP uint32 matrix\n");
		return;
	}
  if(nrhs!=10){
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

  mode_mx = prhs[0];
  image_stack_mx = prhs[1];
  image_stack_scale_mx = prhs[2];
  cat_mx = prhs[3];
  superpixel_2_seg_map_mx = prhs[4];
  pmap_mx = prhs[5];
  mask_label_mx = prhs[6];
  seed_image_mx = prhs[7];
  seed_depth_mx = prhs[8];
  label_new_stack_mx = prhs[9];
 
  mode = (int) * mxGetPr(mode_mx);
	size_0 = mxGetDimensions(image_stack_mx);
	width = size_0[0];
  height = size_0[1];
  depth = size_0[2];
  mexPrintf("%d %d %d\n", width, height, depth);
  //
  // Call the 3D watershed routine
  //
  image_stack = (Stack *) calloc(1, sizeof(Stack));
  image_stack->kind = GREY;
  image_stack->width = width;
  image_stack->height = height;
  image_stack->depth = depth;
  image_stack->array = (uint8 *) mxGetPr(image_stack_mx);

  image_stack_scale = (int) * mxGetPr(image_stack_scale_mx);
  mexPrintf("%d\n", image_stack_scale);
  
  superpixel_maps = (unsigned short int **) calloc(depth, sizeof(unsigned short int *));
  for(i=0; i<depth; i++)
    superpixel_maps[i] = (unsigned short int *) mxGetPr(mxGetCell(cat_mx, i));
	size_1 = mxGetDimensions(mxGetCell(cat_mx, 0));
	original_width = size_1[0];

  superpixel_2_seg_maps = (double **) calloc(depth, sizeof(double *));
  for(i=0; i<depth; i++)
    superpixel_2_seg_maps[i] = mxGetPr(mxGetCell(superpixel_2_seg_map_mx, i));

  pmap = (unsigned int *) mxGetPr(pmap_mx);
  
  mask_label = * (unsigned int *) mxGetPr(mask_label_mx);
  mexPrintf("%u\n", mask_label);
  
  seed_image = (Image *) calloc(1, sizeof(Image));
  seed_image->kind = GREY;
  seed_image->width = width;
  seed_image->height = height;
  seed_image->array = (uint8 *) mxGetPr(seed_image_mx);
  
  seed_depth = * (unsigned int *) mxGetPr(seed_depth_mx);
  
  label_new = (Stack *) calloc(1, sizeof(Stack));
  label_new->kind = GREY;
  label_new->width = width;
  label_new->height = height;
  label_new->depth = depth;
  label_new->array = (uint8 *) mxGetPr(label_new_stack_mx);

  Build_3D_Watershed_Masked_Seeded(image_stack, image_stack_scale, 1, superpixel_maps, original_width, superpixel_2_seg_maps, pmap, mask_label, seed_image, seed_depth, label_new);
  
  free(image_stack);
  free(seed_image);
  free(label_new);
  
  return;
}
