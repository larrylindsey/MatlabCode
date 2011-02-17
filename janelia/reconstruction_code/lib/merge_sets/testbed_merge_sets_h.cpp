
#include <iostream>
#include <ext/hash_map>
#include <vector>

#include <merge_sets_h.h>

using namespace __gnu_cxx;

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

struct equi
{
  bool operator()(const unsigned int s1, const unsigned int s2) const
  {
    return s1==s2;
  }
};

Merge_Sets_H<unsigned int, hash<unsigned int>, equi> m(NULL);

void print_adam(unsigned int n){
  printf("adams:\n");
  unsigned int i;
  for(i=0; i<n; i++){
    printf("%u ", m.adam[i]);
  }
  printf("\n");
}
void print_adam_ids(unsigned int n){
  printf("adam ids:\n");
  unsigned int i;
  for(i=0; i<n; i++){
    printf("%u ", m.get_adam_id(i));
  }
  printf("\n");
}

int main(int argc, char * argv[]){

  unsigned int i;
  for(i=0; i<10; i++)
    m.add_new_set(i);

  print_adam(10);
  print_adam_ids(10);  

  m.merge(1, 3);
  print_adam(10);
  print_adam_ids(10);

  m.merge(0, 2);
  print_adam(10);
  print_adam_ids(10);

  m.merge(3, 2);
  print_adam(10);
  print_adam_ids(10);

  m.merge(3, 2);
  print_adam(10);
  print_adam_ids(10);

  m.update_adams();
  print_adam(10);
  print_adam_ids(10);
  
  return 0;
}
