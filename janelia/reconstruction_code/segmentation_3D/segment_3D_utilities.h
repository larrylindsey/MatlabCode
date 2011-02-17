// Utility functions for 3D segmentation

#ifndef __SEGMENT_3D_UTILITIES__
#define __SEGMENT_3D_UTILITIES__

#include <stack_3D_utilities.h>

extern int stack_width, stack_height, stack_depth, stack_n_pixel, stack_n_pixel_plane;

// After mergers remove boundaries between merged segments.
void remove_merged_boundaries(unsigned int * label_stack, int stack_width, int stack_height, int stack_depth,
                              unsigned int * output_label_stack);

// Draw boundaries between segments if not already present
void draw_segment_boundaries(unsigned int * label_stack, int stack_width, int stack_height, int stack_depth,
                              unsigned int * output_label_stack);

// Check if a pixel is on the boundary of a segment - for MATLAB column major
inline bool check_is_boundary_MATLAB(unsigned int * label_map, unsigned int index){

  bool is_boundary = false;
  unsigned int s = label_map[index];
  unsigned int neighbors[27];
  int n, n_neighbor = get_neighbor_indexes_MATLAB(index, neighbors);
  for(n=0; n<n_neighbor; n++)
    if(label_map[neighbors[n]]!=s){
      is_boundary = true;
      break;
    }
  return is_boundary;
}

// Check if a pixel is on the boundary of a segment
inline bool check_is_boundary(unsigned int * label_map, unsigned int index){

  bool is_boundary = false;
  unsigned int s = label_map[index];
  unsigned int neighbors[27];
  int n, n_neighbor = get_neighbor_indexes(index, neighbors);
  for(n=0; n<n_neighbor; n++)
    if(label_map[neighbors[n]]!=s){
      is_boundary = true;
      break;
    }
  return is_boundary;
}

// Check if a pixel is on the boundary of a segment
inline bool check_is_boundary_6(unsigned int * label_map, unsigned int index){

  bool is_boundary = false;
  unsigned int s = label_map[index];
  unsigned int neighbors[6];
  int n, n_neighbor = get_neighbor_indexes_6(index, neighbors);
  for(n=0; n<n_neighbor; n++)
    if(label_map[neighbors[n]]!=s){
      is_boundary = true;
      break;
    }
  return is_boundary;
}

#endif
