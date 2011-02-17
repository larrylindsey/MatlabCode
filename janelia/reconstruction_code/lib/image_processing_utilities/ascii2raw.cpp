//ascii to raw

#include <iostream>

#include <file_input_output.h>

int main(int argc, char * argv[]){
  void * a;
  iput::io::Type_Element type_element;
  int n_dimension;
  int * dimensions;
  char * file_name;

  { 
    std::cin >> type_element;

    std::cin >> n_dimension;

    int n_element=1;
    dimensions = (int *) new int[n_dimension];
    for(int i=0; i<n_dimension; i++){
      std::cin >> dimensions[i];
      n_element *= dimensions[i];
    }

    // get type of the variable
    unsigned int type_id = iput::io::get_type_id(type_element);
    unsigned int size_element = iput::io::get_size_element(type_element);
    a = (void *) new unsigned char[size_element*n_element];

    if(type_id == iput::io::type_id_char())
      for(int i=0; i<n_element; i++)
	std::cin >>  ((char *)a)[i];

    if(type_id == iput::io::type_id_unsigned_char())
      for(int i=0; i<n_element; i++)
	std::cin >> ((unsigned char *)a)[i];
    
    if(type_id == iput::io::type_id_short_int())
      for(int i=0; i<n_element; i++)
	std::cin >> ((short int *)a)[i];

    if(type_id == iput::io::type_id_unsigned_short_int())
      for(int i=0; i<n_element; i++)
	std::cin >> ((unsigned short int *)a)[i];

    if(type_id == iput::io::type_id_int())
      for(int i=0; i<n_element; i++)
	std::cin >> ((int *)a)[i];

    if(type_id == iput::io::type_id_unsigned_int())
      for(int i=0; i<n_element; i++)
	std::cin >> ((unsigned int *)a)[i];

    if(type_id == iput::io::type_id_long_int())
      for(int i=0; i<n_element; i++)
	std::cin >> (( int *)a)[i];

    if(type_id == iput::io::type_id_unsigned_long_int())
      for(int i=0; i<n_element; i++)
	std::cin >> ((unsigned long int *)a)[i];

    if(type_id == iput::io::type_id_float())
      for(int i=0; i<n_element; i++)
	std::cin >> ((float *)a)[i];

    if(type_id == iput::io::type_id_double())
      for(int i=0; i<n_element; i++)
	std::cin >> ((double *)a)[i];
  }

  file_name = (char *) new char[1];
  file_name[0] = '\0';

  std::cout << iput::io::fwrite_raw_array(file_name, type_element,
					  n_dimension, dimensions, a);
  delete [] file_name;
  delete [] (unsigned char *)a;
}
