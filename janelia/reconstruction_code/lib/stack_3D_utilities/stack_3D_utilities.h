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

#ifndef __STACK_3D_UTILITIES__
#define __STACK_3D_UTILITIES__

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

int stack_width, stack_height, stack_depth, stack_n_pixel, stack_n_pixel_plane;

extern inline int get_neighbor_indexes_MATLAB(unsigned int index, unsigned int * neighbor_index){

  int n_neighbor=0;
  int index_p = index % stack_n_pixel_plane;
  int z = index / stack_n_pixel_plane;
  int x = index_p / stack_height;
  int y = index_p % stack_height;
  if(z>0){
    if(x>0){
      if(y>0)
        neighbor_index[n_neighbor++] = index - stack_n_pixel_plane - stack_height - 1;
      neighbor_index[n_neighbor++] = index - stack_n_pixel_plane - stack_height;
      if(y<stack_height-1)
        neighbor_index[n_neighbor++] = index - stack_n_pixel_plane - stack_height + 1;
    }
    if(y>0)
      neighbor_index[n_neighbor++] = index - stack_n_pixel_plane - 1;
    neighbor_index[n_neighbor++] = index - stack_n_pixel_plane;
    if(y<stack_height-1)
      neighbor_index[n_neighbor++] = index - stack_n_pixel_plane + 1;
    if(x<stack_width-1){
      if(y>0)
        neighbor_index[n_neighbor++] = index - stack_n_pixel_plane + stack_height - 1;
      neighbor_index[n_neighbor++] = index - stack_n_pixel_plane + stack_height;
      if(y<stack_height-1)
        neighbor_index[n_neighbor++] = index - stack_n_pixel_plane + stack_height + 1;
    }
  }
  if(x>0){
    if(y>0)
      neighbor_index[n_neighbor++] = index-stack_height-1;
    neighbor_index[n_neighbor++] = index-stack_height;
    if(y<stack_height-1)
      neighbor_index[n_neighbor++] = index-stack_height+1;
  }
  if(y>0)
    neighbor_index[n_neighbor++] = index-1;
  if(y<stack_height-1)
    neighbor_index[n_neighbor++] = index+1;
  if(x<stack_width-1){
    if(y>0)
      neighbor_index[n_neighbor++] = index+stack_height-1;
    neighbor_index[n_neighbor++] = index+stack_height;
    if(y<stack_height-1)
      neighbor_index[n_neighbor++] = index+stack_height+1;
  }
  if(z<stack_depth-1){
    if(x>0){
      if(y>0)
        neighbor_index[n_neighbor++] = index + stack_n_pixel_plane - stack_height - 1;
      neighbor_index[n_neighbor++] = index + stack_n_pixel_plane - stack_height;
      if(y<stack_height-1)
        neighbor_index[n_neighbor++] = index + stack_n_pixel_plane - stack_height + 1;
    }
    if(y>0)
      neighbor_index[n_neighbor++] = index + stack_n_pixel_plane - 1;
    neighbor_index[n_neighbor++] = index + stack_n_pixel_plane;
    if(y<stack_height-1)
      neighbor_index[n_neighbor++] = index + stack_n_pixel_plane + 1;
    if(x<stack_width-1){
      if(y>0)
        neighbor_index[n_neighbor++] = index + stack_n_pixel_plane + stack_height - 1;
      neighbor_index[n_neighbor++] = index + stack_n_pixel_plane + stack_height;
      if(y<stack_height-1)
        neighbor_index[n_neighbor++] = index + stack_n_pixel_plane + stack_height + 1;
    }
  }
  return n_neighbor;
}


extern inline int get_neighbor_indexes(unsigned int index, unsigned int * neighbor_index){

  int n_neighbor=0;
  int index_p = index % stack_n_pixel_plane;
  int z = index / stack_n_pixel_plane;
  int y = index_p / stack_width;
  int x = index_p % stack_width;
  if(z>0){
    if(y>0){
      if(x>0)
        neighbor_index[n_neighbor++] = index - stack_n_pixel_plane - stack_width - 1;
      neighbor_index[n_neighbor++] = index - stack_n_pixel_plane - stack_width;
      if(x<stack_width-1)
        neighbor_index[n_neighbor++] = index - stack_n_pixel_plane - stack_width + 1;
    }
    if(x>0)
      neighbor_index[n_neighbor++] = index - stack_n_pixel_plane - 1;
    neighbor_index[n_neighbor++] = index - stack_n_pixel_plane;
    if(x<stack_width-1)
      neighbor_index[n_neighbor++] = index - stack_n_pixel_plane + 1;
    if(y<stack_height-1){
      if(x>0)
        neighbor_index[n_neighbor++] = index - stack_n_pixel_plane + stack_width - 1;
      neighbor_index[n_neighbor++] = index - stack_n_pixel_plane + stack_width;
      if(x<stack_width-1)
        neighbor_index[n_neighbor++] = index - stack_n_pixel_plane + stack_width + 1;
    }
  }
  if(y>0){
    if(x>0)
      neighbor_index[n_neighbor++] = index-stack_width-1;
    neighbor_index[n_neighbor++] = index-stack_width;
    if(x<stack_width-1)
      neighbor_index[n_neighbor++] = index-stack_width+1;
  }
  if(x>0)
    neighbor_index[n_neighbor++] = index-1;
  if(x<stack_width-1)
    neighbor_index[n_neighbor++] = index+1;
  if(y<stack_height-1){
    if(x>0)
      neighbor_index[n_neighbor++] = index+stack_width-1;
    neighbor_index[n_neighbor++] = index+stack_width;
    if(x<stack_width-1)
      neighbor_index[n_neighbor++] = index+stack_width+1;
  }
  if(z<stack_depth-1){
    if(y>0){
      if(x>0)
        neighbor_index[n_neighbor++] = index + stack_n_pixel_plane - stack_width - 1;
      neighbor_index[n_neighbor++] = index + stack_n_pixel_plane - stack_width;
      if(x<stack_width-1)
        neighbor_index[n_neighbor++] = index + stack_n_pixel_plane - stack_width + 1;
    }
    if(x>0)
      neighbor_index[n_neighbor++] = index + stack_n_pixel_plane - 1;
    neighbor_index[n_neighbor++] = index + stack_n_pixel_plane;
    if(x<stack_width-1)
      neighbor_index[n_neighbor++] = index + stack_n_pixel_plane + 1;
    if(y<stack_height-1){
      if(x>0)
        neighbor_index[n_neighbor++] = index + stack_n_pixel_plane + stack_width - 1;
      neighbor_index[n_neighbor++] = index + stack_n_pixel_plane + stack_width;
      if(x<stack_width-1)
        neighbor_index[n_neighbor++] = index + stack_n_pixel_plane + stack_width + 1;
    }
  }
  return n_neighbor;
}

extern inline int get_neighbor_indexes_6(unsigned int index, unsigned int * neighbor_index){

  int n_neighbor=0;
  int index_p = index % stack_n_pixel_plane;
  int z = index / stack_n_pixel_plane;
  int y = index_p / stack_width;
  int x = index_p % stack_width;

  if(z>0)
    neighbor_index[n_neighbor++] = index - stack_n_pixel_plane;
  if(z<stack_depth-1)
    neighbor_index[n_neighbor++] = index + stack_n_pixel_plane;
  if(y>0)
    neighbor_index[n_neighbor++] = index - stack_width;
  if(y<stack_height-1)
    neighbor_index[n_neighbor++] = index + stack_width;
  if(x>0)
    neighbor_index[n_neighbor++] = index - 1;
  if(x<stack_width-1)
    neighbor_index[n_neighbor++] = index + 1;
  
  return n_neighbor;
}

// After mergers remove boundaries between merged segments.
void remove_merged_boundaries(unsigned int * label_stack, int stack_width, int stack_height, int stack_depth,
                              unsigned int * output_label_stack);

// Check if a voxel is at the outer shell of the stack
bool check_is_at_outershell(unsigned int index){
  int z = index / stack_n_pixel_plane;
  if(z==0 || z==stack_depth-1)
    return true;
  int index_p = index % stack_n_pixel_plane;
  int x = index_p % stack_width;
  if(x==0 || x==stack_width-1)
    return true;
  int y = index_p / stack_width;
  if(y==0 || y==stack_height-1)
    return true;
  return false;
}

#endif
