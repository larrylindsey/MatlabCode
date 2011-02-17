// Remove merged boundaries
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
//
// v0	02242009	init. code
//

#include <stack_3D_utilities.h>

extern int stack_width, stack_height, stack_depth, stack_n_pixel, stack_n_pixel_plane;

/*
 * input params:
 * 	1. label stack
 * output:
 *	1. label stack with merged boundaries deleted.
 */
void remove_merged_boundaries(unsigned int * label_map, int width, int height, int depth,
                              unsigned int * output_label_map){

  stack_width = width;
  stack_height = height;
  stack_depth = depth;
	stack_n_pixel = stack_height*stack_width*stack_depth;
  stack_n_pixel_plane = stack_width*stack_height;

  {
    unsigned int i, n;
    unsigned char b;
    unsigned int neighbors[27];
    unsigned int *l, *o, s, s1;
    int n_neighbor;
    l = label_map;
    o = output_label_map;
    for(i=0; i<stack_n_pixel; i++, l++,o++){
      if(*l!=0){
        *o=*l;
        continue;
      }

      n_neighbor = get_neighbor_indexes(i, neighbors);
      // find one non-zero label adjacent to the boundary
      s1=0;
      for(n=0; n<n_neighbor; n++){
        s = label_map[neighbors[n]];
        if(s!=0){
          s1=s;
          break;
        }
      }
      // if no non-zero labels then do nothing
      if(s1==0){
        *o=*l;
        continue;
      }
      // if there is one more non-zero label different from s1
      // then retain boundary pixel (=0) otherwise obliterate it.
      *o=s1;
      for(n=0; n<n_neighbor; n++){
        s = label_map[neighbors[n]];
        if(s!=0 && s!=s1){
          *o=0;
          break;
        }
      }
    }
  }

  return;
}
