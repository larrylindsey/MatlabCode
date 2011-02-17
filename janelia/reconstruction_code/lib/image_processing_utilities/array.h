// n-dimensional array
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
//

#ifndef __IPUJ_ARRAY__
#define __IPUJ_ARRAY__

#include <string.h>

#include <file_input_output.h>
#include <ext/hash_map>

namespace iput
{
  template <typename T>
  class Array {
  public:
    Array(){
      buffer = NULL;
      n_dimension = -1;
      n_element = -1;
    }

    Array(const Array <T> & b);
    
    ~Array(){
      if(buffer!=NULL && n_dimension>0)
        delete [] buffer;
    }
    
    T * buffer;
    int n_dimension, n_element;
    std::vector<int> dimensions;

    // input/output functions
    static Array<T> fread_raw_array(char * filename);
    unsigned int fwrite_raw_array(char * filename);

    Array<T> & operator=(const Array<T> & b);
    
  private:
    inline bool is_valid_raw_array_type(io::Type_Element type_element);
    inline io::Type_Element get_type_element(void);
  };
    
  template <typename T>
  Array<T> Array<T>::fread_raw_array(char * filename){
    Array<T> array;

    array.buffer = NULL;
    io::Type_Element type_element;
    int n_dimension, * dimensions;
    void * read_data = io::fread_raw_array(filename, type_element, n_dimension,
                                           dimensions);
    if(iput::io::is_error_type_element(type_element) || read_data==NULL){
      array.n_dimension = -1;
      return array;
    }
    if(!array.is_valid_raw_array_type(type_element)){
      array.n_dimension = -2;
      return array;
    }
    
    array.n_dimension = n_dimension;
    {
      int i;
      array.n_element = 1;
      for(i=0; i<n_dimension; i++){
        array.dimensions.push_back(dimensions[i]);
	array.n_element *= dimensions[i];
      }
    }
    array.buffer = (T *) read_data;
    
    delete [] dimensions;
    return array;
  }

  template <typename T>
  unsigned int Array<T>::fwrite_raw_array(char * filename){
    int * dims = (int *) new int[n_dimension];
    {
      int i;
      for(i=0; i<n_dimension; i++)
        dims[i] = dimensions[i];
    }
    io::Type_Element te = get_type_element();
    unsigned int err =
      iput::io::fwrite_raw_array(filename, te,
                                 n_dimension, dims, (void *) buffer);

    delete [] dims;
    return err;
  }

  template <typename T>
    Array<T> & Array<T>::operator=(const Array<T> & b){
    if(this==&b)
      return *this;

    // check if allocation can be avoided
    bool alloc_needed = false;
    if(n_dimension==b.n_dimension){
      int j;
      std::vector<int>::iterator i;
      for(i=dimensions.begin(), j=0; i!=dimensions.end(); i++, j++)
	if(*i!=b.dimensions[j]){
	  alloc_needed = true;
	  break;
	}
    }
    else
      alloc_needed = true;

    if(alloc_needed){
      if(this->buffer!=NULL)
	delete [] buffer;

      n_dimension = b.n_dimension;
      dimensions = b.dimensions;
      n_element=1;
      for(std::vector<int>::iterator i=dimensions.begin(); i!=dimensions.end();
	  i++)
	n_element*=*i;
      buffer = (T* ) new T[n_element];
    }

    memcpy((void *) buffer, (const void *)(b.buffer), n_element*sizeof(T));

    return *this;
  }

  template <typename T>
    Array<T>::Array(const Array<T> & b){
    this->n_dimension = -1;
    this->n_element = -1;
    this->buffer = NULL;
    *this = b;
  }

  template <>
  inline bool Array<char>::is_valid_raw_array_type(io::Type_Element type_element){
    return io::get_type_id(type_element)==io::type_id_char();
  }
  template <>
  inline bool Array<unsigned char>::is_valid_raw_array_type(io::Type_Element type_element){
    return io::get_type_id(type_element)==io::type_id_unsigned_char();
  }
  template <>
  inline bool Array<short int>::is_valid_raw_array_type(io::Type_Element type_element){
    return io::get_type_id(type_element)==io::type_id_short_int();
  }
  template <>
  inline bool Array<unsigned short int>::is_valid_raw_array_type(io::Type_Element type_element){
    return io::get_type_id(type_element)==io::type_id_unsigned_short_int();
  }
  template <>
  inline bool Array<int>::is_valid_raw_array_type(io::Type_Element type_element){
    return io::get_type_id(type_element)==io::type_id_int();
  }
  template <>
  inline bool Array<unsigned int>::is_valid_raw_array_type(io::Type_Element type_element){
    return io::get_type_id(type_element)==io::type_id_unsigned_int();
  }
  template <>
  inline bool Array<long int>::is_valid_raw_array_type(io::Type_Element type_element){
    return io::get_type_id(type_element)==io::type_id_long_int();
  }
  template <>
  inline bool Array<unsigned long int>::is_valid_raw_array_type(io::Type_Element type_element){
    return io::get_type_id(type_element)==io::type_id_unsigned_long_int();
  }
  template <>
  inline bool Array<float>::is_valid_raw_array_type(io::Type_Element type_element){
    return io::get_type_id(type_element)==io::type_id_float();
  }
  template <>
  inline bool Array<double>::is_valid_raw_array_type(io::Type_Element type_element){
    return io::get_type_id(type_element)==io::type_id_double();
  }

  template <>
  inline io::Type_Element Array<char>::get_type_element(void){
    return io::get_type_element(io::type_id_char(), sizeof(char));
  }
  template <>
  inline io::Type_Element Array<unsigned char>::get_type_element(void){
    return io::get_type_element(io::type_id_unsigned_char(),
                                      sizeof(unsigned char));
  }
  template <>
  inline io::Type_Element Array<short int>::get_type_element(void){
    return io::get_type_element(io::type_id_short_int(),
                                      sizeof(short int));
  }
  template <>
  inline io::Type_Element Array<unsigned short int>::get_type_element(void){
    return io::get_type_element(io::type_id_unsigned_short_int(),
                                      sizeof(unsigned short int));
  }
  template <>
  inline io::Type_Element Array<int>::get_type_element(void){
    return io::get_type_element(io::type_id_int(),
                                      sizeof(int));
  }
  template <>
  inline io::Type_Element Array<unsigned int>::get_type_element(void){
    return io::get_type_element(io::type_id_unsigned_int(),
                                      sizeof(unsigned int));
  }
  template <>
  inline io::Type_Element Array<long int>::get_type_element(void){
    return io::get_type_element(io::type_id_long_int(),
                                      sizeof(long int));
  }
  template <>
  inline io::Type_Element Array<unsigned long int>::get_type_element(void){
    return io::get_type_element(io::type_id_unsigned_long_int(),
                                      sizeof(unsigned long int));
  }
  template <>
  inline io::Type_Element Array<float>::get_type_element(void){
    return io::get_type_element(io::type_id_float(),
                                      sizeof(float));
  }
  template <>
  inline io::Type_Element Array<double>::get_type_element(void){
    return io::get_type_element(io::type_id_double(),
                                      sizeof(double));
  }  
   
}

#endif // __IPUJ_ARRAY__
