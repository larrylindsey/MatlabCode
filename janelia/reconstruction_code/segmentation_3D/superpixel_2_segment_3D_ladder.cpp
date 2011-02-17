// compute watershed hierarchy such that each segment has at least minimum area
// An alternative to Ladder.m designed by Yuriy Mischenko ~ Dec. 2007
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	01042008	init. code
// v1 03112009  stand alone code.
//

#include <math.h>
#include <algorithm>
#include <stdio.h>

#include <image_lib.h>

#include <file_input_output.h>
#include <segment_3D_utilities.h>
#include <stack_3D_utilities.h>
#include <merge_sets.h>

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

extern int stack_width, stack_height, stack_depth, stack_n_pixel, stack_n_pixel_plane;
unsigned int * n_pixel_of_label;
void merge_area(unsigned int s1, unsigned int s2){
  n_pixel_of_label[s1] += n_pixel_of_label[s2];
}

//#define DUMP_LOG

#define ARGC_BOUNDARY_MAP 1
#define ARGC_SUPERPIXEL_MAP 2
#define ARGC_F_THRESHOLD 3
#define ARGC_AREA_THRESHOLD 4
#define ARGC_OUTPUT_LABEL_MAPPING 5
#define ARGC_OUTPUT_LABEL_MAP 6
int main(int argc, char * argv[]){
	if(argc<ARGC_OUTPUT_LABEL_MAPPING){
		printf("Usage superpixel_2_segment_3D_ladder boundary_stack superpixel_stack f_threshold area_threshold superpixel_2_segment_map [output_segment_stack]\n");
		printf("input params:\n");
		printf("\t1. scalar field on which watershed was performed MxNxP uint8 matrix\n");
		printf("\t2. superpixel label map MxNxP uint32 matrix\n");
		printf("\t3. threshold on the scalar field below which all watershed segments must be merged uint8\n");
		printf("\t4. threshold on the area - minimum number pixels in a segment for it to exist on its own. If violated, the segment is merged uint32.\n");
		printf("\t5. mapping from input segment labels to merged segment labels\n");
		printf("\t6. merged watershed label map MxNxP uint32 (optional)\n");
		return 1;
	}

  printf("\nSTART: superpixel_2_segment_3D_ladder_b\n");
  Stack *stack= Read_Stack(argv[ARGC_BOUNDARY_MAP]);
	printf("height:%d width:%d depth:%d\n", stack->height, stack->width, stack->depth);

  #ifdef DUMP_LOG
  FILE * fout = fopen("ladder.log", "wt");
  #endif
  
	stack_width = stack->width;
  stack_height = stack->height;
  stack_depth = stack->depth;
	stack_n_pixel = stack_height*stack_width*stack_depth;
  stack_n_pixel_plane = stack_width*stack_height;
  printf("stack_width=%d, stack_height=%d, stack_depth=%d, stack_n_pixel=%d, stack_n_pixel_plane=%d\n",
            stack_width, stack_height, stack_depth, stack_n_pixel, stack_n_pixel_plane);

	unsigned int * superpixel_label_map;
  {
    iput::io::Type_Element type_element;
    int n_dimension, * dimensions;
    void * read_data = iput::io::fread_raw_array(argv[ARGC_SUPERPIXEL_MAP],
                                                 type_element, n_dimension,
                                                 dimensions);
    if(iput::io::is_error_type_element(type_element) || read_data==NULL){
      printf("Could not open superpixel map\n");
      return -1;
    }
    if(n_dimension!=3){
      printf("Not a 3D segment stack\n");
      return -1;
    }
    if(dimensions[0]!=stack_width || dimensions[1]!=stack_height || dimensions[2]!=stack_depth){
      printf("Segment map's dimensions don't match with those of image stack\n");
      printf("Segment map's dimensions: width:%d height:%d depth:%d\n",
             dimensions[0], dimensions[1], dimensions[2]);
      return -1;
    }
    superpixel_label_map = (unsigned int *) read_data;
  }
 
	unsigned char * boundary_value = (unsigned char *) stack->array;
	unsigned char f_threshold = (unsigned char) atoi(argv[ARGC_F_THRESHOLD]);
  printf("f_threshold = %d\n", f_threshold);
	unsigned int area_threshold = (unsigned int) atoi(argv[ARGC_AREA_THRESHOLD]);
  printf("area_threshold = %d\n", area_threshold);
  
  // Sort the boundary pixels in increasing order of value.
  // Sort by indexing (linear time).
  // Find boundary pixels of a particular value and form lists.
  printf("Collecting statistics of boundary values and superpixels ... ");
  unsigned int n_boundary_pixel_value[256];
  unsigned int n_boundary_pixel_value_cumsum[256];
  unsigned int max_superpixel_id = 0;
  unsigned int n_boundary_pixel = 0;
  {
    unsigned int i;
    for(i=0; i<256; i++)
      n_boundary_pixel_value[i]=0;
    unsigned int neighbors[27];
    int n_neighbor, n;
    unsigned int * s = superpixel_label_map;
    unsigned char * b = boundary_value;
    bool is_boundary;
    for(i=0; i<stack_n_pixel; i++, s++, b++){
      max_superpixel_id = MAX(max_superpixel_id, *s);

      is_boundary = check_is_boundary(superpixel_label_map, i);
      if(is_boundary){
        n_boundary_pixel_value[*b]++;
        n_boundary_pixel++;
      }
    }
    int current_sum = 0;
    for(i=0; i<256; i++){
      current_sum = n_boundary_pixel_value_cumsum[i] = n_boundary_pixel_value[i]+current_sum;
      n_boundary_pixel_value_cumsum[i]--;
    }

    printf("n_boundary_pixel = %d\n", n_boundary_pixel);
    printf("max_superpixel_id = %d\n", max_superpixel_id);
/*    
    // DEBUG
    for(i=0; i<256; i++)
      printf("%d ", n_boundary_pixel_value_cumsum[i]);
    printf("\n");
*/  
  }
  
  unsigned int * buffer_1 = (unsigned int *) new unsigned int[max_superpixel_id+1];
	n_pixel_of_label = buffer_1;
  {
    int i;
    unsigned int * l=n_pixel_of_label;
    for(i=0; i<=max_superpixel_id; i++, l++){
      *l=0;
    }
	}
  printf("done.\n");
  
  printf("Sorting the boundary values and collecting other stats ... ");
  unsigned int * buffer_2 = (unsigned int *) new unsigned int[n_boundary_pixel];
  unsigned int * boundary_value_index = buffer_2;
  {
    unsigned int i;
    unsigned char * b = boundary_value;
    unsigned int * s = superpixel_label_map;
    bool is_boundary;
    for(i=0; i<stack_n_pixel; i++, b++, s++){
      is_boundary = check_is_boundary_6(superpixel_label_map, i);
      if(is_boundary){
        boundary_value_index[n_boundary_pixel_value_cumsum[*b]] = i;
        n_boundary_pixel_value_cumsum[*b]--;
      }
      else
        n_pixel_of_label[*s]++;
    }

    /*
    // DEBUG
    for(i=0; i<stack_n_pixel; i++)
      printf("%d ", boundary_value_index[i]);
    printf("\n");
    */

/*    
    // DEBUG
    for(i=0; i<max_superpixel_id+1; i++)
      printf("%d %d\n", i, n_pixel_of_label[i]);
    printf("\n");
*/  
  }
  printf("done.\n");
  
  {
    /*
    // DEBUG
    unsigned int neighbor_index[27];
    int n_neighbor;
    int i;
    n_neighbor = get_neighbor_indexes(13, neighbor_index);
    printf("%d\n", n_neighbor);
    //return;
    for(i=0; i<n_neighbor; i++)
      printf("%d ", neighbor_index[i]);
    printf("\n");
    */
  }
  
  // Loop through the boundary values making mergers.
  // Merge all superpixels with boundary less than threshold or if
  // one of them has area less than threshold.
  Merge_Sets superpixel_sets(max_superpixel_id, merge_area);
  
  {
    int i, n1, n2;
    unsigned char b;
    unsigned int b_i, neighbors[27];
    unsigned int s1, s2, s3;
    unsigned int adam_s1, adam_s2, area_1, area_2, temp_adams[27];
    int n_neighbor;
    for(i=0; i<n_boundary_pixel; i++){
      if(i%100000==0)
        printf("*", i);
      if(i%4000000==0)
        printf("\n");
      fflush(stdout);
      
      b_i = boundary_value_index[i];
      b = boundary_value[b_i];

      if(b<f_threshold){
        #ifdef DUMP_LOG
        fprintf(fout, "boundary merge: %d\n", b);
        #endif
        n_neighbor = get_neighbor_indexes_6(b_i, neighbors);
        for(n1=0; n1<n_neighbor; n1++){
          s1 = superpixel_label_map[neighbors[n1]];
          if(s1==0)
            continue;
          for(n2=0; n2<n_neighbor; n2++){
            if(n1==n2)
              continue;
            s2 = superpixel_label_map[neighbors[n2]];
            if(s2==0)
              continue;
            if(s1==s2)
              continue;
            if(s1>s2){
              s3=s1;
              s1=s2;
              s2=s3;
            }
            superpixel_sets.merge(s1,s2);
          }
        }
      }
      else{
        n_neighbor = get_neighbor_indexes(b_i, neighbors);
        for(n1=0; n1<n_neighbor; n1++){
          s1 = superpixel_label_map[neighbors[n1]];
          if(s1==0)
            continue;
          adam_s1 = superpixel_sets.get_adam(s1);
          area_1 = n_pixel_of_label[adam_s1];
          for(n2=0; n2<n1; n2++){
            s2 = superpixel_label_map[neighbors[n2]];
            if(s2==0)
              continue;
            adam_s2 = superpixel_sets.get_adam(s2);
            area_2 = n_pixel_of_label[adam_s2];
            if(adam_s1==adam_s2)
              continue;
            if(area_1<area_threshold ||
               area_2<area_threshold){
              // one of the segment areas is too small - merge.
              #ifdef DUMP_LOG
              fprintf(fout, "before area merge: s1:%u ads1:%u s2:%u ads2:%u as1:%u as2:%u\n",
                      s1, adam_s1, s2, adam_s2, area_1, area_2);
              #endif
              if(s1>s2)
                superpixel_sets.merge(s2,s1);
              else
                superpixel_sets.merge(s1,s2);
              
              adam_s1 = superpixel_sets.get_adam(s1);
              area_1 = n_pixel_of_label[adam_s1];
              #ifdef DUMP_LOG
              fprintf(fout, "after area merge: s1:%u ads1:%u s2:%u ads2:%u as1:%u as2:%u\n",
                      s1, adam_s1, s2, adam_s2, area_1, area_2);
              #endif
            }
          }
        }
      }
    }
  }

  printf("\nFinalizing the label mappings ... ");
  superpixel_sets.update_adams();
  printf("done\n");

  {
    int size_element = sizeof(unsigned int);
    iput::io::Type_Element type_element = iput::io::get_type_element(
      iput::io::type_id_unsigned_int(), size_element);
    int n_dimension = 1;
    int dimensions[1] = {max_superpixel_id+1};
    int err;
    {
      unsigned int * adam_c = (unsigned int *) new unsigned int[superpixel_sets.adam.size()];
      int i;
      for(i=0; i<superpixel_sets.adam.size(); i++)
        adam_c[i] = superpixel_sets.adam[i];
      err = iput::io::fwrite_raw_array(argv[ARGC_OUTPUT_LABEL_MAPPING],
                                       type_element, n_dimension, dimensions,
                                       (void *) adam_c);
      delete [] adam_c;
    }
    if(err!=0){
      printf("ERROR superpixel_2_segment_3D_ladder_b: could not write label mapping [%d]\n", err);
      delete [] buffer_1;
      delete [] buffer_2;
      return -2;
    }
  }

  if(argc>ARGC_OUTPUT_LABEL_MAP){
    unsigned int * segment_label_map = (unsigned int *)
      new unsigned int[stack_n_pixel];
    int i;
    unsigned int * l = segment_label_map;
    unsigned int * s = superpixel_label_map;
    for(i=0; i<stack_n_pixel; i++, s++, l++)
      *l = superpixel_sets.adam[*s];

    int size_element = sizeof(unsigned int);
    iput::io::Type_Element type_element = iput::io::get_type_element(
      iput::io::type_id_unsigned_int(), size_element);
    int n_dimension = 3;
    int dimensions[3] = {stack_width, stack_height, stack_depth};
    int err;
    err = iput::io::fwrite_raw_array(argv[ARGC_OUTPUT_LABEL_MAP],
                                     type_element, n_dimension, dimensions,
                                     (void *) segment_label_map);
    if(err!=0){
      printf("ERROR superpixel_2_segment_3D_ladder_b: could not write label map [%d]\n", err);
      delete [] segment_label_map;
      delete [] buffer_1;
      delete [] buffer_2;
      return -2;
    }
    
    delete [] segment_label_map;
  }
  printf("\nSTOP: superpixel_2_segment_3D_ladder_b\n");
  Kill_Stack(stack);
  delete [] buffer_1;
  delete [] buffer_2;
  return 0;
}
