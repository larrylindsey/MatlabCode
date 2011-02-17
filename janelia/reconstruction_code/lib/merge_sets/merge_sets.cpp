/*
 * Class for merge sets defined over a superset of N elements. 
 */

#include <iostream>
#include <merge_sets.h>

Merge_Sets::Merge_Sets(unsigned int n=0, void (*f)(unsigned int, unsigned int)=NULL){
  // this allows referring to the elements directly without first
  // subtracting 1
  n_element = n+1;
  
  adam.reserve(n_element);
  for(int i=0; i<n_element; i++)
    adam.push_back(i);
  call_on_merge=f;
}

Merge_Sets::Merge_Sets(unsigned int n=0){
  Merge_Sets(n, NULL);
}

Merge_Sets::~Merge_Sets(){
  return;
}

// Add new sets - (may be) useful for incremental algorithms
void Merge_Sets::add_new_sets(unsigned int n_new_element){
  adam.reserve(adam.size()+n_new_element);
  for(int i=0; i<n_new_element; i++){
    adam.push_back(adam.size());
  }

  n_element += n_new_element;
}

void Merge_Sets::print_adams(void){
  update_adams();
  std::vector<unsigned int>::iterator i;
  int j;
/*  for(i=adam.begin(), j=0; i!=adam.end(); i++, j++){
    std::cout << *i << " ";
    if((j+1)%20==0)
      std::cout << "\n";
  }
  std::cout << "\n";
 */
}

// Merge the sets containing elements e1 and e2
void Merge_Sets::merge(unsigned int e1, unsigned int e2){
  unsigned int e;

  unsigned int a1 = adam[e1];
  while(adam[a1]!=a1)
    a1=adam[a1];
  while(e1!=a1){
    e = adam[e1];
    adam[e1]=a1;
    e1 = e;
  }

  while(e2!=adam[e2]){
    e = adam[e2];
    adam[e2]=a1;
    e2 = e;
  }
  unsigned int a2 = adam[e2];
  adam[e2]=a1;

  if(call_on_merge!=NULL)
    call_on_merge(a1, a2);
}

// Get adam of element e - performs an update to make future queries fast
unsigned int Merge_Sets::get_adam(unsigned int e){
  update_adams(e);
  return adam[e];
}

// Update adams of element e
void Merge_Sets::update_adams(unsigned int e){
  unsigned a = adam[e], e1;
  if(a==e || a==adam[a])
    return; // e is either an adam or pointing to an adam - nothing to do.
  // get the adam
  while(adam[a]!=a)
    a=adam[a];
  // set all elements in the chain from e to a to point to a.
  while(e!=a){
    e1 = adam[e];
    adam[e]=a;
    e = e1;
  }
}

// Update adams of all elements
void Merge_Sets::update_adams(void){
  for(int i=0; i<n_element; i++){
    unsigned a = adam[i], e, e1;
    if(a==i || a==adam[a])
      continue; // i is either an adam or pointing to an adam - nothing to do.
    // get the adam
    while(a!=adam[a])
      a=adam[a];
    // set all elements in the chain from e to a to point to a.
    e = i;
    while(e!=a){
      e1 = adam[e];
      adam[e]=a;
      e = e1;
    }
  }
}
