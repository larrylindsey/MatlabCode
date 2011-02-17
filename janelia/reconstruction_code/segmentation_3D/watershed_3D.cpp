/*
 * Thin wrapper for 3D watershed implementation by Eugene Myers, JFRC, HHMI.
 * License: Contact Eugene Myers.
 *
 * Shiv N. Vitaladevuni,
 * Janelia Farm Research Campus, HHMI.
 */
#include <stdio.h>
#include <stdlib.h>

#include <image_lib.h>
#include <utilities.h>
#include <water_shed.h>

#include <file_input_output.h>
#include <filter.h>

int main(int argc, char * argv[]){
  static char * Spec [] = { "[!V] [-v <int>]",
			    " [-{f|filter} <string>]", 
			    " <boundary_map:string>", 
			    " <output_segment_map:string>",
			    NULL  };

  Process_Arguments(argc, argv, Spec, 0);
  
  int is_verbose = 0;
  if(Is_Arg_Matched("-V"))
    is_verbose = 255;
  if(Is_Arg_Matched("-v @"))
    is_verbose = Get_Int_Arg("-v @");
  
  if(is_verbose>0)
    printf("START: %s\n", argv[0]);

  Stack *stack= Read_Stack(Get_String_Arg("boundary_map"));
  iput::Array<unsigned char> boundary_map;
  boundary_map.n_dimension = 3;
  boundary_map.dimensions.push_back(stack->width);
  boundary_map.dimensions.push_back(stack->height);
  boundary_map.dimensions.push_back(stack->depth);
  boundary_map.n_element = stack->width*stack->height*stack->depth;
  boundary_map.buffer = (unsigned char *) stack->array;
  int nvoxel = boundary_map.n_element;
  if(is_verbose)
    printf("height:%d width:%d depth:%d\n", 
	   stack->height, stack->width, stack->depth);

  if(Is_Arg_Matched("-f @")){
    char * filter_param = Get_String_Arg("-f @");
    if(is_verbose>1)
      std::cout << "Applying filter " << filter_param << '\n';
    iput::Array<unsigned char> b;
    iput::filter(boundary_map, b, filter_param, is_verbose-1);
    boundary_map = b;
  }
  stack->array = (unsigned char *) boundary_map.buffer;

  Stack * watershed_label;
  watershed_label = Build_3D_Watershed3(stack, 1);
  Kill_Stack(stack);

  // Save as unsigned int array
  if(sizeof(int)!=sizeof(unsigned int)){
    printf("Assumption that sizeof(int)==sizeof(unsigned int) failed!\n");
    return -1;
  }
  {
    int * buffer = (int *) watershed_label->array;
    for(int i=0; i<boundary_map.n_element; i++)
      if(buffer[i]<0)
	buffer[i]=0;
  }
  unsigned int size_element = sizeof(unsigned int);
  iput::io::Type_Element type_element = 
    iput::io::get_type_element(iput::io::type_id_unsigned_int(),
                               size_element);
  int n_dimension = 3;
  int dimensions[3];
  dimensions[0] = stack->width;
  dimensions[1] = stack->height;
  dimensions[2] = stack->depth;
  int result = iput::io::fwrite_raw_array(Get_String_Arg("output_segment_map"),
                                          type_element, n_dimension,
                                          dimensions,
                                          (void *) watershed_label->array);
  Kill_Stack(watershed_label);
  boundary_map.buffer = NULL; // freeing via stack variable
  
  if(is_verbose>0)
    printf("STOP: %s\n", argv[0]);
  
  return 0;
}
