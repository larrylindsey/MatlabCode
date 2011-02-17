/*
 * Thin wraper for 3D watershed implementation by Eugene Myers, JFRC, HHMI.
 * License: Contact Eugene Myers.
 *
 * Shiv N. Vitaladevuni
 * Janelia Farm Research Campus,
 * Howard Hughes Medical Institute.
 */ 
#include <stdio.h>
#include <stdlib.h>
#include <image_lib.h>
#include <water_shed.h>

#define ARGC_FORMAT 1
#define ARGC_STACK 2
#define ARGC_MASK 3
int main(int argc, char * argv[]){
  /* Initialize stack */
  printf("Reading stack .. ");
  Stack *stack= Read_Stack(argv[ARGC_STACK]);
  printf("done.\n");
  printf("height=%d, width=%d, depth=%d\n", stack->height, stack->width, stack->depth);
  
  int nvoxel = stack->height * stack->width * stack->depth;
  int i;
  uint8 *array = (uint8 *)stack->array;
  for (i = 0; i < nvoxel; i++){
//    array[i] = 0xFF - array[i];
    array[i] = array[i]<90?0:array[i];
  }
  
  int format = atoi(argv[ARGC_FORMAT]);
  int out_file_argv = 0;
  Watershed_3D * watershed;
  switch(format){
    case ARGC_STACK:
    default:
      printf("Performing vanilla watershed .. ");
      watershed = Build_3D_Watershed(stack, 1);
      out_file_argv = 3;
      break;
    case ARGC_MASK:
      printf("Performing seeded watershed .. ");
      Stack * mask = Read_Stack(argv[2]);
      watershed = Build_3D_Watershed_Masked(stack, 1, mask, 1);
      uint8 * mask_array = (uint8 *) mask->array;
      for (i = 0; i < nvoxel; i++)
        array[i] = (mask_array[i]>0)*(0xFF - array[i]);
      Kill_Stack(mask);
      out_file_argv = 4;
      break;
  }
  printf("done\n");
  
  printf("Dumping watershed result stack .. ");
  Stack * watershed_stack = Color_Watersheds_3D(watershed, stack);
  Write_Stack(argv[out_file_argv], watershed_stack);
  printf("done\n");
  
  Kill_Stack(watershed_stack);
  Kill_Stack(stack);
  
  return 0;
}
