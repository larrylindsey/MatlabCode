// Generate a blend of segmentation label color map and the grayscale stack for
// display purposes.
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	03122009	init. code
//

#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <image_lib.h>

#include <segment_3D_utilities.h>
#include <stack_3D_utilities.h>
#include <merge_sets.h>
#include <file_input_output.h>
#include <colormaps.h>

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

#define ARGC_STACK_NAME_FORMAT 1
#define ARGC_SUPERPIXEL_NAME_FORMAT 2
#define ARGC_SEGMENT_NAME_FORMAT 3
#define ARGC_OUTPUT_BLEND_TIF 4
#define ARGC_OUTPUT_SCALE 5
#define ARGC_FIRST_SECTION 6
int main(int argc, char * argv[]){
	if(argc<ARGC_OUTPUT_BLEND_TIF){
		printf("Usage blend_stitch_substack_segment substack_name_format superpixel_name_format segment_name_format output_blend_tif section_sequence\n");
		printf("input params:\n");
		printf("\t1. Format for the image sub-stack file names.\n");
		printf("\t2. Format for the superpixel file names.\n");
		printf("\t3. Format for segment file names.\n");
		printf("\t4. blended image-segment stack RGB uint8\n");
    printf("\t5. scale of output image in x-y plane.\n");
    printf("\t6. Sequence of begin and end sections of the sub-stacks.\n");
		return 1;
	}

  printf("\nSTART: blend_segment_stack_b\n");
  printf("stack file name format: %s\n",
         argv[ARGC_STACK_NAME_FORMAT]);
  printf("superpixel file name format: %s\n",
         argv[ARGC_SUPERPIXEL_NAME_FORMAT]);
  printf("segment file name format: %s\n",
         argv[ARGC_SEGMENT_NAME_FORMAT]);
  int n_substack = 0;
  {
    printf("Sub-stack begin and end:\n");
    int i;
    for(i=ARGC_FIRST_SECTION; i<argc; i+=2){
      printf("%d %d\n", atoi(argv[i]), atoi(argv[i+1]));
      n_substack++;
    }
  }

  int output_scale = atoi(argv[ARGC_OUTPUT_SCALE]);
  
  // Segment labels are hashed using multiplication method with bit XOR.
  // Let m be a bitwise-random mask. Each bit has a probability p of
  // being 1 and (1-p) of being 0.
  // Let A be a "nice" irrational number, e.g., (sqrt(5)-1)/2
  // h(l) = floor(n_color * ((l XOR m)*A mod 1))
  unsigned int bitwise_random_mask;
  {
    int i;
    // generate bitwise-random mask
    int p = 50;
    srand(time(NULL));
    bitwise_random_mask = (rand()%100)<p;
    printf("Bitwise-random mask: %d", bitwise_random_mask);
    int m_1;
    for(i=0; i<32;i++){
      m_1 = (rand()%100)<p;
      printf("%d", m_1);
      bitwise_random_mask <<= m_1;
    }
    printf("\n");
  }
  double A = 0.5*(sqrt(5)-1);

  Stack * blend_stack = NULL;
  unsigned int stack_width_s;
  unsigned int stack_height_s;
  unsigned int stack_n_pixel_plane_s;

  char * file_name = (char *) new char
    [MAX(MAX(strlen(argv[ARGC_SUPERPIXEL_NAME_FORMAT]),
             strlen(argv[ARGC_SEGMENT_NAME_FORMAT])),
         strlen(argv[ARGC_STACK_NAME_FORMAT]))+50];
  int substack_id;
  for(substack_id=0; substack_id<n_substack; substack_id++){
    int argc_substack = substack_id*2 + ARGC_FIRST_SECTION;
    int start_section = atoi(argv[argc_substack]);
    int end_section = atoi(argv[argc_substack+1]);
  
    sprintf(file_name, argv[ARGC_STACK_NAME_FORMAT],
            start_section, end_section);
    printf("Reading substack %s\n", file_name);
    Stack *stack= Read_Stack(file_name);
    int stack_width = stack->width;
    int stack_height = stack->height;
    int stack_depth = stack->depth;
    int stack_n_pixel = stack_height*stack_width*stack_depth;
    int stack_n_pixel_plane = stack_width*stack_height;
    printf("stack_width=%d, stack_height=%d, stack_depth=%d, stack_n_pixel=%d, stack_n_pixel_plane=%d, total_depth=%d\n",
           stack_width, stack_height, stack_depth, stack_n_pixel, stack_n_pixel_plane, atoi(argv[argc-1]) - atoi(argv[ARGC_FIRST_SECTION])+1 );

    if(blend_stack==NULL){
      stack_width_s = stack_width/output_scale;
      stack_height_s = stack_height/output_scale;
      stack_n_pixel_plane_s = stack_width_s*stack_height_s;
      blend_stack = Make_Stack(COLOR, stack_width_s,
                               stack_height_s,
                               atoi(argv[argc-1]) -
                               atoi(argv[ARGC_FIRST_SECTION])+1);
    }
    
    printf("Reading superpixel for substack %d\n", substack_id);
    unsigned int * superpixel_map;
    {
      sprintf(file_name, argv[ARGC_SUPERPIXEL_NAME_FORMAT],
              start_section, end_section);
      iput::io::Type_Element type_element;
      int n_dimension, * dimensions;
      void * read_data = iput::io::fread_raw_array(file_name,
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
      if(dimensions[0]!=stack_width || dimensions[1]!=stack_height
         || dimensions[2]!=stack_depth){
        printf("\nERROR: Segment map's dimensions don't match with those of image stack\n");
        printf("\nERROR: Segment map's dimensions: width:%d height:%d depth:%d\n",
               dimensions[0], dimensions[1], dimensions[2]);
        return -1;
      }
      superpixel_map = (unsigned int *) read_data;
      delete [] dimensions;
    }

    unsigned int max_superpixel_label = 0;
    {
      int i;
      unsigned int * l = superpixel_map;
      for(i=0; i<stack_n_pixel; i++, l++)
        max_superpixel_label = MAX(max_superpixel_label, *l);
    }

    printf("Reading segment mapping for substack %d\n", substack_id);
    unsigned int * label_mapping;
    {
      sprintf(file_name, argv[ARGC_SEGMENT_NAME_FORMAT],
              start_section, end_section);
      iput::io::Type_Element type_element;
      int n_dimension, * dimensions;
      void * read_data = iput::io::fread_raw_array(file_name,
                                                   type_element, n_dimension,
                                                   dimensions);
      if(iput::io::is_error_type_element(type_element) || read_data==NULL){
        printf("\nERROR: Could not open label mapping\n");
        return -1;
      }
      if(n_dimension!=1){
        printf("\nERROR: Not a 1D label mapping\n");
        return -1;
      }
      if(dimensions[0]!=max_superpixel_label+1){
        printf("\nERROR: Label mapping's dimensions %d don't match with max segment label %d\n", dimensions[0], max_superpixel_label+1);
        return -1;
      }
      label_mapping = (unsigned int *) read_data;
      delete [] dimensions;
    }

    printf("Generating blended substack %d ...\n", substack_id);
    unsigned char * b = (unsigned char *) blend_stack->array +
      3 * stack_n_pixel_plane_s *
      (start_section - atoi(argv[ARGC_FIRST_SECTION]));
    {
      int d, x_s, y_s, x_o, y_o, h;
      unsigned char * c;
      // original and scaled indexes into stack
      unsigned char * gray_stack = ((unsigned char *) stack->array);
      unsigned int index_o, index_s;
      for(d=0; d<stack_depth; d++){
        for(y_o=0, y_s=0; y_s<stack_height_s;
            y_s++, y_o+=output_scale){
          for(x_o=0, x_s=0; x_s<stack_width_s;
              x_s++, x_o+=output_scale){
            index_o = d*stack_n_pixel_plane + y_o*stack_width + x_o;
            index_s = d*stack_n_pixel_plane_s + y_s*(stack_width_s)
              + x_s;
            h = floor(fmod(A*(double)(label_mapping[superpixel_map[index_o]]
                                      ^ bitwise_random_mask), 1)
                      * iput::display::colormap_26.n_color);
            c = iput::display::colormap_26.map + 3*h;
            b[index_s*3] = ((255-(int)(gray_stack[index_o])) * (*c))/255;
            b[index_s*3+1] = ((255-(int)(gray_stack[index_o])) * (*(c+1)))/255;
            b[index_s*3+2] = ((255-(int)(gray_stack[index_o])) * (*(c+2)))/255;
          }
        }
      }
    }
    Kill_Stack(stack);
    delete [] superpixel_map;
    delete [] label_mapping;
    printf("done.\n");
  }

  printf("Dumping tif ... ");
  fflush(stdout);
  Write_Stack(argv[ARGC_OUTPUT_BLEND_TIF], blend_stack);
  printf("done.\n");

  delete [] file_name;
  printf("\nSTOP: blend_segment_stack_b\n");
  return 0;
}
