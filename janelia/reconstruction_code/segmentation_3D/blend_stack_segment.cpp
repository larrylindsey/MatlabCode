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

#include <image_lib.h>

#include <segment_3D_utilities.h>
#include <stack_3D_utilities.h>
#include <merge_sets.h>
#include <file_input_output.h>
#include <colormaps.h>

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

int main(int argc, char * argv[]){
  int ARGC_IMAGE_STACK = 1;
  int ARGC_SEGMENT_MAP = 2;
  int ARGC_SEGMENT_MAPPING = 3;
  int ARGC_OUTPUT_BLEND_TIF = 4;

  printf("no. of arguments: %d\n", argc);

	if(argc<=ARGC_SEGMENT_MAPPING){
		printf("Usage blend_stack_segment image_stack segment_stack segment_mapping output_blend_tif\n");
		printf("input params:\n");
		printf("\t1. scalar field on which segmentation was performed MxNxP uint8 matrix\n");
		printf("\t2. segment label map MxNxP uint32 matrix\n");
		printf("\t3. segment mapping: map labels in segment_stack to merged labels Qx1 uint32\n");
		printf("\t4. blended image-segment stack RGB uint8\n");
		return 1;
	}

  if(argc<=ARGC_OUTPUT_BLEND_TIF){
    ARGC_OUTPUT_BLEND_TIF = ARGC_SEGMENT_MAPPING;
    ARGC_SEGMENT_MAPPING = -1;
  }

  printf("\nSTART: blend_segment_stack_b\n");
  Stack *stack= Read_Stack(argv[ARGC_IMAGE_STACK]);
	printf("height:%d width:%d depth:%d\n", stack->height, stack->width, stack->depth);

	int stack_width = stack->width;
  int stack_height = stack->height;
  int stack_depth = stack->depth;
	int stack_n_pixel = stack_height*stack_width*stack_depth;
  int stack_n_pixel_plane = stack_width*stack_height;
  printf("stack_width=%d, stack_height=%d, stack_depth=%d, stack_n_pixel=%d, stack_n_pixel_plane=%d\n",
            stack_width, stack_height, stack_depth, stack_n_pixel, stack_n_pixel_plane);

	unsigned int * segment_map;
  {
    iput::io::Type_Element type_element;
    int n_dimension, * dimensions;
    void * read_data = iput::io::fread_raw_array(argv[ARGC_SEGMENT_MAP],
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
    if(dimensions[0]!=stack_width || dimensions[1]!=stack_height ||
       dimensions[2]!=stack_depth){
      printf("\nERROR: Segment map's dimensions don't match with those of image stack\n");
      printf("\nERROR: Segment map's dimensions: width:%d height:%d depth:%d\n",
             dimensions[0], dimensions[1], dimensions[2]);
      return -1;
    }
    segment_map = (unsigned int *) read_data;
    delete [] dimensions;
  }

  unsigned int max_segment_label = 0;
  {
    int i;
    unsigned int * l = segment_map;
    for(i=0; i<stack_n_pixel; i++, l++)
      max_segment_label = MAX(max_segment_label, *l);
  }

  unsigned int * label_mapping = NULL;
  if(ARGC_SEGMENT_MAPPING>0){
    iput::io::Type_Element type_element;
    int n_dimension, * dimensions;
    void * read_data = iput::io::fread_raw_array(argv[ARGC_SEGMENT_MAPPING],
                                                 type_element, n_dimension,
                                                 dimensions);
    if(iput::io::is_error_type_element(type_element) || read_data==NULL){
      printf("\nERROR: Could not open label mapping\n");
      return -1;
    }
    if(n_dimension>2){
      printf("\nERROR: Not a 1D label mapping\n");
      return -1;
    }
    if(dimensions[0]!=max_segment_label+1){
      printf("\nERROR: Label mapping's dimensions don't match with max segment label\n");
      return -1;
    }
    label_mapping = (unsigned int *) read_data;
    delete [] dimensions;
  }
  
  // Generate blended stack
  // Segment labels are hashed using multiplication method with bit XOR.
  // Let m be a bitwise-random mask. Each bit has a probability p of being 1 and
  // (1-p) of being 0.
  // Let A be a "nice" irrational number, e.g., (sqrt(5)-1)/2
  // h(l) = floor(n_color * ((l XOR m)*A mod 1))
  Stack * blend_stack = Make_Stack(COLOR, stack_width, stack_height, stack_depth);
  {
    unsigned char * g = (unsigned char *) stack->array;
    unsigned char * b = (unsigned char *) blend_stack->array;
    unsigned int m; // bitwise-random mask
    int i;
    // generate bitwise-random mask
    int p = 50;
    srand(time(NULL));
    m = (rand()%100)<p;
    printf("Bitwise-random mask: %d", m);
    int m_1;
    for(i=0; i<32;i++){
      m_1 = (rand()%100)<p;
      printf("%d", m_1);
      m <<= m_1;
    }
    printf("\n");

    double A = 0.5*(sqrt(5)-1);
    
    int h;
    unsigned char * c;
    unsigned int * l = segment_map;
    unsigned int s;
    for(i=0; i<stack_n_pixel; i++, l++, b+=3, g++){
      if(label_mapping!=NULL)
        s = label_mapping[*l];
      else
        s = *l;
      h = floor(fmod(A*(double)(s ^ m), 1) *
                iput::display::colormap_26.n_color);
      c = iput::display::colormap_26.map + 3*h;
      *b = ((255-(int)(*g)) * (*c))/255;
      *(b+1) = ((255-(int)(*g)) * (*(c+1)))/255;
      *(b+2) = ((255-(int)(*g)) * (*(c+2)))/255;
    }
  }
  Write_Stack(argv[ARGC_OUTPUT_BLEND_TIF], blend_stack);

  printf("\nSTOP: blend_segment_stack_b\n");
  Kill_Stack(stack);
  delete [] segment_map;
  if(label_mapping!=NULL)
    delete [] label_mapping;
  return 0;
}
