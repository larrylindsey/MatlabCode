// utilities for file input and output of image processing datastructures
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI.
//
// v0  03122009  init code
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "file_input_output.h"

namespace iput
{
namespace io
{
  
  const Type_Element rv_no_error = 0;
  const Type_Element rv_invalid_parameter = 0x80000001;
  const Type_Element rv_could_not_open = 0x80000002;
  const Type_Element rv_could_not_write = 0x80000003;
  const Type_Element rv_could_not_read = 0x80000004;
  const Type_Element rv_could_not_allocate = 0x80000005;

  Type_Element type_id_char(){ return 0x00000000;}
  Type_Element type_id_unsigned_char(){ return 0x00000100;}
  Type_Element type_id_short_int(){ return 0x00000200;}
  Type_Element type_id_unsigned_short_int(){ return 0x00000300;}
  Type_Element type_id_int(){return 0x00000400;}
  Type_Element type_id_unsigned_int(){return 0x00000500;}
  Type_Element type_id_long_int(){return 0x00000600;}
  Type_Element type_id_unsigned_long_int(){return 0x00000700;}
  Type_Element type_id_float(){return 0x00000800;}
  Type_Element type_id_double(){return 0x00000900;}

  unsigned int get_type_id(Type_Element type_element){
    return type_element & 0xffffff00;
  }
  
  unsigned int get_size_element(Type_Element type_element){
    return type_element & 0x000000ff;
  }
  
  Type_Element get_type_element(unsigned int type_id,
                                unsigned int size_element){
    return (Type_Element) (type_id | size_element);
  }
  
  bool is_error_type_element(Type_Element type_element){
    return (type_element & 0x80000000) > 0;
  }
  
  unsigned int get_error_type_element(Type_Element type_element){
    return type_element & 0x000000ff;
  }
  
  /*
    type_element = [data_type 3 bytes][#bytes per element 1byte]
    data_type:
    0 char
    1 unsigned char
    2 short int
    3 unsigned short int
    4 int
    5 unsigned int
    6 long int
    7 unsigned long int
    8 float
    9 double
   */
  Type_Element fwrite_raw_array(char * file_name, Type_Element type_element,
                       int n_dimension, int * dimensions, void * data){
    int d;
    if(n_dimension<=0)
      return rv_invalid_parameter;
    for(d=0; d<n_dimension; d++)
      if(dimensions[d]<0)
        return rv_invalid_parameter;
    
    FILE * fout;
    if(strlen(file_name)==0)
      fout = stdout;
    else
      fout = fopen(file_name, "wb");

    if(fout==NULL)
      return rv_could_not_open;
    
    if(fwrite(&type_element, sizeof(unsigned int), 1, fout)!=1){
      fclose(fout);
      return rv_could_not_write;
    }
    if(fwrite(&n_dimension, sizeof(int), 1, fout)!=1){
      fclose(fout);
      return rv_could_not_write;
    }
    if(fwrite(dimensions, sizeof(int), n_dimension, fout)!=n_dimension){
      fclose(fout);
      return rv_could_not_write;
    }
    int n_element = dimensions[0];
    for(d=1; d<n_dimension; d++)
      n_element *= dimensions[d];
    unsigned int size_element = get_size_element(type_element);
    if(fwrite(data, size_element, n_element, fout)!=n_element){
      fclose(fout);
      return rv_could_not_write;
    }
    fclose(fout);
    return rv_no_error;
  }

  void * fread_raw_array(char * file_name, Type_Element & type_element,
                         int & n_dimension,int * & dimensions){
    FILE * fin;
    if(strlen(file_name)==0)
      fin = stdin;
    else
      fin = fopen(file_name, "rb");

    if(fin==NULL){
      type_element = rv_could_not_open;
      return NULL;
    }
  
    if(fread(&type_element, sizeof(unsigned int), 1, fin)!=1){
      type_element = rv_could_not_read;
      fclose(fin);
      return NULL;
    }
    if(fread(&n_dimension, sizeof(int), 1, fin)!=1){
      type_element = rv_could_not_read;
      fclose(fin);
      return NULL;
    }
    dimensions = (int *) new int[n_dimension];
    if(dimensions==NULL){
      type_element = -1;
      fclose(fin);
      return NULL;
    }
    if(fread(dimensions, sizeof(int), n_dimension, fin)!=n_dimension){
      type_element = rv_could_not_read;
      fclose(fin);
      return NULL;
    }
  
    int d;
    int n_element = dimensions[0];
    for(d=1; d<n_dimension; d++)
      n_element *= dimensions[d];
    if(n_element<=0){
      return NULL;
    }
    unsigned int size_element = get_size_element(type_element);
    char * data = (char *) new char*[size_element*n_element];
    if(data==NULL){
      type_element = rv_could_not_allocate;
      return NULL;
    }
    if(fread(data, size_element, n_element, fin)!=n_element){
      type_element = rv_could_not_read;
      delete [] data;
      fclose(fin);
      return NULL;
    }
    fclose(fin);
    return (void *) data;
  }

  int check_file_input_output(char * file_name, int n_dimension, 
			      int * dimensions, bool is_random_data){
    {
#define DEALLOC {delete [] data; if(data_read!=NULL) delete [] data_read;}
      int i;
      unsigned int size_element = sizeof(double);
      unsigned int type_element = 0x00000900 | size_element;
      int n_element = dimensions[0];
      for(i=1; i<n_dimension; i++)
        n_element *= dimensions[i];
      double * data;
      if(n_element==0)
        data = NULL;
      else
        data = (double *) new double[n_element];
      if(n_element!=0 && data==NULL){
        printf("ERROR test_file_input_output: allocating testing buffer\n");
        return rv_could_not_allocate;
      }
      printf("Generating random data ...");
      for(i=0; i<n_element; i++)
	if(is_random_data)
	  data[i] = (double) (rand() % 1000);
	else
	  data[i] = i;
      printf("done.\n");
      
      unsigned int err;
      if((err=fwrite_raw_array(file_name, type_element, n_dimension,
                               dimensions, (void *)data))!=0){
        printf("ERROR test_file_input_output: writing test file [%d]\n", err);
        delete [] data;
        return err;
      }

      unsigned int type_element_read;
      int n_dimension_read;
      int * dimensions_read;
      double * data_read = NULL;
      data_read = (double *) fread_raw_array(file_name, type_element_read,
                                             n_dimension_read, dimensions_read);
      if(data!=NULL && data_read==NULL){
        printf("ERROR test_file_input_output: in reading test file: [%d]\n",
               type_element_read);
        DEALLOC;
        return type_element_read;
      }

      if(type_element!=type_element_read){
        printf("ERROR test_file_input_output: type_element != type_element_read\n");
        DEALLOC;
        return -1;
      }
      if(n_dimension!=n_dimension_read){
        printf("ERROR test_file_input_output: n_dimension != n_dimension_read\n");
        DEALLOC;
        return -1;
      }
      for(i=0; i<n_dimension; i++){
        if(dimensions[i]!=dimensions_read[i]){
          printf("ERROR test_file_input_output: dimensions[%d] != dimensions_read[%d]\n", i, i);
          DEALLOC;
          return -1;
        }
      }
      for(i=0; i<n_element; i++){
        if(data[i]!=data_read[i]){
          printf("ERROR test_file_input_output: data[%d] != data_read[%d]\n", i, i);
          DEALLOC;
          return -1;
        }
      }
    }

    return 0;
  }

  int test_file_input_output(char * file_name){
    int net_result = 0;
    {
      // Test one dimensional array
      int n_dimension = 1;
      int dimensions[1] = {986};
      int result;
      result = check_file_input_output(file_name, n_dimension, dimensions, true);
      if(result==0)
        printf("Result for test on %dx1 double array successful.\n", dimensions[0]);
      else
        printf("ERROR: Result for test on %dx1 double array unsuccessful.\n", dimensions[0]);
      net_result = result!=0?result:net_result;
    }
    {
      // Test two dimensional array
      int n_dimension = 2;
      int dimensions[2] = {432, 391};
      int result;
      result = check_file_input_output(file_name, n_dimension, dimensions, true);
      if(result==0)
        printf("Result for test on %dx%d double array successful.\n", dimensions[1],
               dimensions[0]);
      else
        printf("ERROR: Result for test on %dx%d double array unsuccessful.\n", dimensions[1],
               dimensions[0]);
      net_result = result!=0?result:net_result;
    }      
    {
      // Test three dimensional array
      int n_dimension = 3;
      int dimensions[3] = {432, 391, 10};
      int result;
      result = check_file_input_output(file_name, n_dimension, dimensions, true);
      if(result==0)
        printf("Result for test on %dx%dx%d double array successful.\n", dimensions[1],
               dimensions[0], dimensions[2]);
      else
        printf("ERROR: Result for test on %dx%dx%d double array unsuccessful.\n", dimensions[1],
               dimensions[0], dimensions[2]);
      net_result = result!=0?result:net_result;
    }
    {
      // Test three dimensional array
      int n_dimension = 3;
      int dimensions[3] = {2, 3, 2};
      int result;
      result = check_file_input_output(file_name, n_dimension, dimensions, false);
      if(result==0)
        printf("Result for test on %dx%dx%d double array successful.\n", dimensions[1],
               dimensions[0], dimensions[2]);
      else
        printf("ERROR: Result for test on %dx%dx%d double array unsuccessful.\n", dimensions[1],
               dimensions[0], dimensions[2]);
      net_result = result!=0?result:net_result;
    }
    return net_result;
  }
}
} // end of namespace iput

