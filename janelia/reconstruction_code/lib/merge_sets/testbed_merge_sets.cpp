
#include <iostream>
#include <vector>

#include <merge_sets.h>

std::vector<int> sizes;
static int n_element;
static void merge_size(unsigned int a1, unsigned int a2){
  sizes[a1] += sizes[a2];
}
static void print_sizes(){
  for(int i=0; i<sizes.size(); i++){
    std::cout << sizes[i] << " ";
    if((i+1)%20==0)
      std::cout << "\n";
  }
  std::cout << "\n";
}
  
int main(int argc, char * argv[]){

  n_element = 10;

  for(int i=0; i<=n_element; i++)
    sizes.push_back(1);
  
  Merge_Sets s(n_element, merge_size);
  s.print_adams();
  print_sizes();
  
  std::cout << "s.merge(0,1);\n";
  s.merge(0,1);
  s.print_adams();
  print_sizes();
  
  std::cout << "s.merge(2,3);\n";
  s.merge(2,3);
  s.print_adams();
  print_sizes();
  
  std::cout << "s.merge(1,2);\n";
  s.merge(1,2);
  s.print_adams();
  print_sizes();

  std::cout << "s.add_new_sets(10);\n";
  for(int i=0; i<10; i++)
    sizes.push_back(1);
  s.add_new_sets(10);
  s.print_adams();
  print_sizes();
  
  return 0;
}
