/*
 * Class for merge sets defined over a superset of N elements. 
 */

#ifndef __MERGE_SETS_H__
#define __MERGE_SETS_H__

#include <stdio.h>
#include <ext/hash_map>

using namespace __gnu_cxx;

/*
  Example of template:
  T: unsigned int
  H: hash<unsigned int>
  E: equi where
  struct equi{ bool operator()(const unsigned int s1, const unsigned int s2) const
  {return s1==s2;}};
*/

template <class T, class H, class E>
class Merge_Sets_H{
public:
  // Create superset.
  // Also call a function call_on_merge when merging sets
  // with adams a1 and a2
  Merge_Sets_H(void (*f) (T, T)){ call_on_merge=f; n_element=0;}

  // Destructor
  ~Merge_Sets_H();

  // Add new set
  void add_new_set(T new_element);

  // Add new set if it doesn't exist already - (may be) useful for
  // incremental algorithms
  void add_new_set_inc(T new_element);

  // Merge sets containing elements e1 and e2.
  void merge(T e1, T e2);

  // Get adam for element e - performs update for faster queries in future
  T get_adam(T e);
  unsigned int get_adam_id(T e);
  
  // Update adams for all or specific element e
  void update_adams();
  void update_adams(T e);
  
  // Private data
  unsigned int n_element;
  typedef hash_map<T, T, H, E> Adam_Hash;
  Adam_Hash adam;
  typedef typename Adam_Hash::iterator Adam_Hash_Iterator;
  hash_map<T, unsigned int, H, E> adam_id;
  hash_map<T, unsigned int, H, E> n_member;
  void (*call_on_merge)(T, T);
};

template <class T, class H, class E>
Merge_Sets_H<T,H,E>::~Merge_Sets_H(){
  return;
}

template <class T, class H, class E>
void Merge_Sets_H<T,H,E>::add_new_set(T new_element){
  adam[new_element] = new_element;
  n_element++;
  adam_id[new_element] = n_element;
  n_member[new_element] = 1;
  return;
}

template <class T, class H, class E>
void Merge_Sets_H<T,H,E>::add_new_set_inc(T new_element){
  if(adam[new_element]!=0)
    return;
  add_new_set(new_element);
  return;
}

template <class T, class H, class E>
void Merge_Sets_H<T,H,E>::merge(T e1, T e2){
  T a1 = e1;
  while(a1!=adam[a1])
    a1 = adam[a1];
  unsigned int i1 = adam_id[a1];
  
  T a2 = e2;
  while(a2!=adam[a2])
    a2 = adam[a2];
  unsigned int i2 = adam_id[a2];
  
  if(n_member[a1]<n_member[a2]){
    T b = e1;
    T c;
    while(b!=a1){
      c = adam[b];
      adam[b] = a1;
      adam_id[b] = i1;
      b = c;
    }
    
    b = e2;
    while(b!=a2){
      c = adam[b];
      adam[b] = a1;
      adam_id[b] = i1;
      b = c;
    }
    adam[b] = a1;
    adam_id[b] = i1;

    n_member[a1] += n_member[a2];
    return;
  }

  T b = e2;
  T c;
  while(b!=a2){
    c = adam[b];
    adam[b] = a2;
    adam_id[b] = i2;
    b = c;
  }
    
  b = e1;
  while(b!=a1){
    c = adam[b];
    adam[b] = a2;
    adam_id[b] = i2;
    b = c;
  }
  adam[b] = a2;
  adam_id[b] = i2;

  n_member[a2] += n_member[a1];
  return;
}

template <class T, class H, class E>
T Merge_Sets_H<T,H,E>::get_adam(T e){
  T a = e;
  while(a!=adam[a])
    a = adam[a];
  T b;
  while(e!=a){
    b = adam[e];
    adam[e] = a;
    e = b;
  }
  return a;
}

template <class T, class H, class E>
unsigned int Merge_Sets_H<T,H,E>::get_adam_id(T e){
  return adam_id[get_adam(e)];
}

template <class T, class H, class E>
void Merge_Sets_H<T,H,E>::update_adams(T e){
  T a = e;
  while(a!=adam[a])
    a = adam[a];
  T b;
  while(e!=a){
    b = adam[e];
    adam[e] = a;
    e = b;
  }
}

template <class T, class H, class E>
void Merge_Sets_H<T,H,E>::update_adams(void){
  Adam_Hash_Iterator it;
  T a, b, e;
  for(it=adam.begin(); it!=adam.end(); it++){
    e = (*it).first;
    a = e;
    while(a!=adam[a])
      a = adam[a];
    while(e!=a){
      b = adam[e];
      adam[e] = a;
      e = b;
    }
  }
}
#endif
