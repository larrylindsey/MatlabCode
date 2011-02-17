// utilities for file input and output of image processing datastructures
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI.
//
// v0  03122009  init code
//

#ifndef __IPUJ_FILE_INPUT_OUTPUT__
#define __IPUJ_FILE_INPUT_OUTPUT__

namespace iput
{
namespace io  
{
  /*
    type_element = [data_type 3 bytes][#bytes per element 1byte]
    data_type:
  */

  typedef unsigned int Type_Element;
  Type_Element type_id_char();
  Type_Element type_id_unsigned_char();
  Type_Element type_id_short_int();
  Type_Element type_id_unsigned_short_int();
  Type_Element type_id_int();
  Type_Element type_id_unsigned_int();
  Type_Element type_id_long_int();
  Type_Element type_id_unsigned_long_int();
  Type_Element type_id_float();
  Type_Element type_id_double();


  unsigned int get_type_id(Type_Element type_element);
  unsigned int get_size_element(Type_Element type_element);
  Type_Element get_type_element(unsigned int type_id,
                                unsigned int size_element);

  bool is_error_type_element(Type_Element type_element);
  unsigned int get_error_type_element(Type_Element type_element);
  
  // Write a n_dimension array of elements of size_element bytes each to a file
  // named file_name from a buffer pointed to by data.
  unsigned int fwrite_raw_array(char * file_name, Type_Element type_element,
                                int n_dimension, int * dimensions,
                                void * data);

  // Read a n_dimension array of elements of size_element bytes each from file
  // named file_name. Memory is allocated and pointed to by data.
  void * fread_raw_array(char * file_name, Type_Element & type_element,
                         int & n_dimension, int * & dimensions);

  // Test the input output routines with a temporary file named file_name.
  int test_file_input_output(char * file_name);
}
}
#endif //__IPU_FILE_INPUT_OUTPUT__
