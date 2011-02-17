// get segment boundary sets for a given segmentation and boundary map. This can be
// be used for training segmentation algorithms.
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	08272008	init. code
//

#include <mex.h>
#include <math.h>
#include <algorithm>

#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

#define MAX_3_3(p,H) MAX(MAX(MAX(MAX(*(p-1-H),*(p-1)),MAX(*(p-1+H),*(p+H))),MAX(MAX(*(p+1+H),*(p+1)),MAX(*(p+1-H),*(p-H)))),*(p))
#define MAX_2_2(p,H) MAX(MAX(*(p),*(p+1)),MAX(*(p+H),*(p+H+1)))

#define NUM_ORIENT_BINS		180

class Boundary{
  public:
    int label_0, label_1;
    int x,y;
    double f;
};
bool operator<(const Boundary& a, const Boundary& b){
  return (a.label_0<b.label_0)||((a.label_0==b.label_0)&&(a.label_1<b.label_1))||
  ((a.label_0==b.label_0)&&(a.label_1==b.label_1)&&(a.x<b.x))||
  ((a.label_0==b.label_0)&&(a.label_1==b.label_1)&&(a.x==b.x)&&(a.y<b.y))||
  ((a.label_0==b.label_0)&&(a.label_1==b.label_1)&&(a.x==b.x)&&(a.y==b.y)&&(a.f<b.f));
}


/*
    mexPrintf("Usage seg_boundary_sets = get_segment_boundary_sets(segment_label_map, boundary_map);\n");
    mexPrintf("input params:\n");
    mexPrintf("\t1. segment label map MxN matrix\n");
    mexPrintf("\t2. scalar field of boundary values MxN matrix\n");
    mexPrintf("\t3. scalar field of boundary mask {0,1} MxN matrix\n");
    mexPrintf("output:\n");
    mexPrintf("\t1. structure array of boundary element sets. One set for each pair of\n");
    mexPrintf("\t   adjacent segments. Structure fields: segment_label_0, segment_label_1\n");
    mexPrintf("\t   and boundary_set - a [3xn_i] matrix [x, y, b(x,y)].\n");
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){
  
  if(nrhs==0){
    if(nlhs==1){
      plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
      *mxGetPr(plhs[0]) = 1;
      return;
    }
    mexPrintf("Usage seg_boundary_sets = get_segment_boundary_sets(segment_label_map, boundary_map);\n");
    mexPrintf("input params:\n");
    mexPrintf("\t1. segment label map MxN matrix\n");
    mexPrintf("\t2. scalar field of boundary values MxN matrix\n");
    mexPrintf("\t3. scalar field of boundary mask {0,1} MxN matrix\n");
    mexPrintf("output:\n");
    mexPrintf("\t1. structure array of boundary element sets. One set for each pair of\n");
    mexPrintf("\t   adjacent segments. Structure fields: segment_label_0, segment_label_1\n");
    mexPrintf("\t   and boundary_set - a [3xn_i] matrix [x, y, b(x,y)].\n");
    return;
  }
  if(nrhs!=3){
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }
  if(nlhs>1){
    mexErrMsgTxt("Wrong number of outputs\n");
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
  int n_pixel = height*width;
  double * seg_label_map = mxGetPr(prhs[0]);
  
  double * watershed_field = mxGetPr(prhs[1]);
  double * boundary_mask = mxGetPr(prhs[2]);
  
  int * n_pixel_of_label = (int *) new int[n_pixel];
  int * l=n_pixel_of_label;
  for(int i=0; i<n_pixel; i++, l++)
    *l=0;
  
  int n_boundary_elements=0;
  double * ws_l = seg_label_map;
  double * b_m = boundary_mask;
  for(int i=0; i<n_pixel; i++, ws_l++, b_m++){
    if((int)*ws_l==0 && (*b_m)==1){
      n_boundary_elements++;
    }
  }
  
  const char * field_names[] = {"segment_label_0", "segment_label_1", "boundary_set"};
  if(n_boundary_elements==0){
    const int temp [] = {1};
    plhs[0] = mxCreateStructArray(1, temp, 3, field_names);
    return;
  }

  /////
  // Construct a list of boundary elements
  /////
  // Padding to take care of image edges
  double * padded_watershed_label_map = (double *) new double[(height+2)*(width+2)];
  double * p_ws_l = padded_watershed_label_map;
  for(int i=0; i<(height+2)*(width+2); i++, p_ws_l++)
    *p_ws_l=0;
  
  p_ws_l = padded_watershed_label_map+height+2+1;
  ws_l = seg_label_map;
  for(int x=0; x<width; x++, p_ws_l+=2){
    for(int y=0; y<height; y++, ws_l++,p_ws_l++){
      *p_ws_l=*ws_l;
    }
  }
  
  // Collect elements
  double * f_map = watershed_field;
  Boundary * boundary_list = (Boundary *) new Boundary[n_boundary_elements*4];
  p_ws_l = padded_watershed_label_map+height+2+1;
  b_m = boundary_mask;
  int j=0;
  int height_padded = height+2;
  int label_0, label_1;
  for(int x=0; x<width; x++, p_ws_l+=2){
    for(int y=0; y<height; y++, f_map++, p_ws_l++, b_m++){
      if((*p_ws_l)==0 && (*b_m)==1){
        label_0=(int)*(p_ws_l-1); label_1=(int)*(p_ws_l+1);
        if(label_0<label_1){
          boundary_list[j].label_0=label_0; 			boundary_list[j].label_1=label_1;
        }
        else{
          boundary_list[j].label_0=label_1; 			boundary_list[j].label_1=label_0;
        }
        boundary_list[j].f=*f_map;
        boundary_list[j].x = x;
        boundary_list[j].y = y;
        j++;
        
        label_0=(int)*(p_ws_l-height_padded); label_1=(int)*(p_ws_l+height_padded);
        if(label_0<label_1){
          boundary_list[j].label_0=label_0; 			boundary_list[j].label_1=label_1;
        }
        else{
          boundary_list[j].label_0=label_1; 			boundary_list[j].label_1=label_0;
        }
        boundary_list[j].f=*f_map;
        boundary_list[j].x = x;
        boundary_list[j].y = y;
        j++;
        
        label_0=(int)*(p_ws_l-height_padded-1); label_1=(int)*(p_ws_l+height_padded+1);
        if(label_0<label_1){
          boundary_list[j].label_0=label_0; 			boundary_list[j].label_1=label_1;
        }
        else{
          boundary_list[j].label_0=label_1; 			boundary_list[j].label_1=label_0;
        }
        boundary_list[j].f=*f_map;
        boundary_list[j].x = x;
        boundary_list[j].y = y;
        j++;
        
        label_0=(int)*(p_ws_l-height_padded+1); label_1=(int)*(p_ws_l+height_padded-1);
        if(label_0<label_1){
          boundary_list[j].label_0=label_0; 			boundary_list[j].label_1=label_1;
        }
        else{
          boundary_list[j].label_0=label_1; 			boundary_list[j].label_1=label_0;
        }
        boundary_list[j].f=*f_map;
        boundary_list[j].x = x;
        boundary_list[j].y = y;
        j++;
      }
    }
  }
  
  // Get the list of boundary elements to be involved in the watershed
  std::sort(boundary_list, boundary_list+n_boundary_elements*4-1);
  
  int * boundary_lengths = (int *) new int[n_boundary_elements];
  int n_adj_seg_pair = 0;
  int length_boundary_list_sorted=0;
  Boundary * b_l = boundary_list;
  int prev_label_0=0, prev_label_1=0;
  int prev_x=-1, prev_y=-1;
  int boundary_length=0;
  for(int i=0; i<n_boundary_elements*4; i++, b_l++){
    if(b_l->label_0<=0 || b_l->label_1<=0)
      continue;
    
    if(b_l->label_1!=prev_label_1 || b_l->label_0!=prev_label_0){
      if(boundary_length!=0){
        boundary_lengths[n_adj_seg_pair] = boundary_length;
        n_adj_seg_pair++;
      }
      
      prev_label_0 = b_l->label_0;
      prev_label_1 = b_l->label_1;
      prev_x = b_l->x;
      prev_y = b_l->y;
      boundary_length = 1;
    }
    else{
      if(b_l->y!=prev_y || b_l->x!=prev_x){
        boundary_length++;
        prev_x = b_l->x;
        prev_y = b_l->y;
      }
    }
  }
  boundary_lengths[n_adj_seg_pair] = boundary_length;
  n_adj_seg_pair++;
  
  int temp_size[1];
  temp_size[0] = n_adj_seg_pair;
  mxArray * seg_boundary_sets = mxCreateStructArray(1, temp_size, 3, field_names);
  for(int i=0; i<n_adj_seg_pair; i++){
    mxSetFieldByNumber(seg_boundary_sets, i, 0, mxCreateDoubleMatrix(1, 1, mxREAL));
    mxSetFieldByNumber(seg_boundary_sets, i, 1, mxCreateDoubleMatrix(1, 1, mxREAL));
    mxSetFieldByNumber(seg_boundary_sets, i, 2, mxCreateDoubleMatrix(3, boundary_lengths[i], mxREAL));
  }
  
  b_l = boundary_list;
  prev_label_0=0, prev_label_1=0;
  prev_x=-1, prev_y=-1;
  boundary_length=0;
  int pair_id = 0;
  double * boundary_set = NULL;
  mxArray * boundary_set_mx;
  for(int i=0; i<n_boundary_elements*4; i++, b_l++){
    if(b_l->label_0<=0 || b_l->label_1<=0)
      continue;
    
    if(b_l->label_1!=prev_label_1 || b_l->label_0!=prev_label_0){
      if(boundary_length!=0){
        pair_id++;
      }
      
      prev_label_0 = b_l->label_0;
      prev_label_1 = b_l->label_1;
      prev_x = b_l->x;
      prev_y = b_l->y;
      //mexPrintf("%d %d\n", n_adj_seg_pair, pair_id);
      boundary_set_mx = mxGetFieldByNumber(seg_boundary_sets, pair_id, 2);
      if(boundary_set_mx==NULL){
        mexPrintf("Boundary set pointer is null!!\n");
        return;
      }
      else{
        boundary_set = mxGetPr(boundary_set_mx);
      }
      boundary_set[0] = b_l->x;
      boundary_set[1] = b_l->y;
      boundary_set[2] = b_l->f;
      *(mxGetPr(mxGetFieldByNumber(seg_boundary_sets, pair_id, 0))) = b_l->label_0;
      *(mxGetPr(mxGetFieldByNumber(seg_boundary_sets, pair_id, 1))) = b_l->label_1;
      boundary_length = 1;
    }
    else{
      if(b_l->y!=prev_y || b_l->x!=prev_x){
        boundary_set[3*boundary_length  ] = b_l->x;
        boundary_set[3*boundary_length+1] = b_l->y;
        boundary_set[3*boundary_length+2] = b_l->f;
        boundary_length++;
        prev_x = b_l->x;
        prev_y = b_l->y;
      }
    }
  }
  
  plhs[0] = seg_boundary_sets;
  
  delete [] n_pixel_of_label;
  delete [] padded_watershed_label_map;
  delete [] boundary_list;
  delete [] boundary_lengths;
  
  return;
}


