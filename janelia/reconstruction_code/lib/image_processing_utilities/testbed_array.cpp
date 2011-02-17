// testbed for Array class

#include <array.h>
#include <iostream>
#include <stdio.h>

int main(int argc, char * argv[]){

  {
    iput::Array<int> a;
    a.n_dimension = 1;
    a.dimensions.push_back(5);
    a.buffer = (int *) new int[5];
    {
      int i;
      for(i=0; i<5; i++)
        a.buffer[i] = i+100;
    }

    char filename[1000] =
      "/groups/chklovskii/home/vitaladevunis/temp/array_test.raw";
    printf("Write status: %u\n", a.fwrite_raw_array(filename));

    iput::Array<int> b = iput::Array<int>::fread_raw_array(filename);
    std::cout << "n_dimension: " << b.n_dimension << '\n';
    std::cout << "dimensions: ";
    copy(b.dimensions.begin(), b.dimensions.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << '\n';
    {
      int n_element = b.dimensions[0];
      for(int i=1; i<b.n_dimension; i++)
        n_element *= b.dimensions[i];
      for(int i=0; i<n_element; i++)
        std::cout << b.buffer[i] << ' ';
      std::cout << '\n';
    }
  }

  {
    iput::Array<unsigned int> a;
    a.n_dimension = 2;
    a.dimensions.push_back(5);
    a.dimensions.push_back(3);
    a.buffer = (unsigned int *) new unsigned int[5*3];
    {
      int i, j;
      for(i=0; i<5; i++)
        for(j=0; j<3; j++)
        a.buffer[j*5+i] = j*1000+i;
    }

    char filename[1000] =
      "/groups/chklovskii/home/vitaladevunis/temp/array_test.raw";
    printf("Write status: %u\n", a.fwrite_raw_array(filename));

    iput::Array<unsigned int> b = iput::Array<unsigned int>::fread_raw_array(filename);
    std::cout << "n_dimension: " << b.n_dimension << '\n';
    std::cout << "dimensions: ";
    copy(b.dimensions.begin(), b.dimensions.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << '\n';
    {
      int n_element = b.dimensions[0];
      for(int i=1; i<b.n_dimension; i++)
        n_element *= b.dimensions[i];
      for(int i=0; i<n_element; i++)
        std::cout << b.buffer[i] << ' ';
      std::cout << '\n';
    }
  }
  
  return 0;
}
