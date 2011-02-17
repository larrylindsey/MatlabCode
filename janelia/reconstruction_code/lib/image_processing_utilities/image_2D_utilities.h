// General utility routines for 3D stacks
// (1) Get 26-connected neighbors' indexes.
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// Matlab version called compute_segmentation_hierarchy_from_watershed_with_min_area.m
//
// v0	01042008	init. code
//

#ifndef __IMAGE_2D_UTILITIES__
#define __IMAGE_2D_UTILITIES__

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

namespace iput
{
  int image_width, image_height, image_n_pixel;

#ifdef __MATLAB_MEX__
  extern inline int get_neighbor_indexes(unsigned int index, unsigned int * neighbor_index){

    int n_neighbor=0;
    int x = index / image_height;
    int y = index % image_height;

    if(x>0){
      if(y>0)
        neighbor_index[n_neighbor++] = index-image_height-1;
      neighbor_index[n_neighbor++] = index-image_height;
      if(y<image_height-1)
        neighbor_index[n_neighbor++] = index-image_height+1;
    }
    if(y>0)
      neighbor_index[n_neighbor++] = index-1;
    if(y<image_height-1)
      neighbor_index[n_neighbor++] = index+1;
    if(x<image_width-1){
      if(y>0)
        neighbor_index[n_neighbor++] = index+image_height-1;
      neighbor_index[n_neighbor++] = index+image_height;
      if(y<image_height-1)
        neighbor_index[n_neighbor++] = index+image_height+1;
    }

    return n_neighbor;
  }
#endif

#ifndef __MATLAB_MEX__
  extern inline int get_neighbor_indexes(unsigned int index, unsigned int * neighbor_index){

    int n_neighbor=0;
    int y = index / image_width;
    int x = index % image_width;

    if(y>0){
      if(x>0)
        neighbor_index[n_neighbor++] = index-image_width-1;
      neighbor_index[n_neighbor++] = index-image_width;
      if(x<image_width-1)
        neighbor_index[n_neighbor++] = index-image_width+1;
    }
    if(x>0)
      neighbor_index[n_neighbor++] = index-1;
    if(x<image_width-1)
      neighbor_index[n_neighbor++] = index+1;
    if(y<image_height-1){
      if(x>0)
        neighbor_index[n_neighbor++] = index+image_width-1;
      neighbor_index[n_neighbor++] = index+image_width;
      if(x<image_width-1)
        neighbor_index[n_neighbor++] = index+image_width+1;
    }

    return n_neighbor;
  }
#endif
  
// After mergers remove boundaries between merged segments.
  void remove_merged_boundaries_2D(
    unsigned int * label_stack,
    int image_width, int image_height, int image_depth,
    unsigned int * output_label_stack);
}
#endif
