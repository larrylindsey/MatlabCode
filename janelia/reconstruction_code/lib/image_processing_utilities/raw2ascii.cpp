//raw file to ascii

#include <iostream>

#include <file_input_output.h>


int main(int argc, char * argv[]){
  void * a;
  iput::io::Type_Element type_element;
  int n_dimension;
  int * dimensions;
  char * file_name;

  file_name = (char *) new char[1];
  file_name[0] = '\0';

  a = iput::io::fread_raw_array(file_name, type_element,
				n_dimension, dimensions);

  std::cout << type_element << ' ' <<  n_dimension << ' ';
  int n_element;
  {
    n_element=1;
    for(int i=0; i<n_dimension; i++){
      std::cout << dimensions[i] << ' ';
      n_element*=dimensions[i];
    }
  }

  { 
    // get type of the variable
    unsigned int type_id = iput::io::get_type_id(type_element);

    if(type_id == iput::io::type_id_char())
      for(int i=0; i<n_element; i++)
	std::cout << ((char *)a)[i] << ' ';

    if(type_id == iput::io::type_id_unsigned_char())
      for(int i=0; i<n_element; i++)
	std::cout << ((unsigned char *)a)[i] << ' ';
    
    if(type_id == iput::io::type_id_short_int())
      for(int i=0; i<n_element; i++)
	std::cout << ((short int *)a)[i] << ' ';

    if(type_id == iput::io::type_id_unsigned_short_int())
      for(int i=0; i<n_element; i++)
	std::cout << ((unsigned short int *)a)[i] << ' ';

    if(type_id == iput::io::type_id_int())
      for(int i=0; i<n_element; i++)
	std::cout << ((int *)a)[i] << ' ';

    if(type_id == iput::io::type_id_unsigned_int())
      for(int i=0; i<n_element; i++)
	std::cout << ((unsigned int *)a)[i] << ' ';

    if(type_id == iput::io::type_id_long_int())
      for(int i=0; i<n_element; i++)
	std::cout << (( int *)a)[i] << ' ';

    if(type_id == iput::io::type_id_unsigned_long_int())
      for(int i=0; i<n_element; i++)
	std::cout << ((unsigned long int *)a)[i] << ' ';

    if(type_id == iput::io::type_id_float())
      for(int i=0; i<n_element; i++)
	std::cout << ((float *)a)[i] << ' ';

    if(type_id == iput::io::type_id_double())
      for(int i=0; i<n_element; i++)
	std::cout << ((double *)a)[i] << ' ';
    
    std::cout << "\n";
  }

  delete [] file_name;
  delete [] (unsigned char *) a;
}
