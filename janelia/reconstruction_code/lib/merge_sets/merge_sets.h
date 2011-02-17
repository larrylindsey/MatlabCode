/*
 * Class for merge sets defined over a superset of N elements. 
 */

#ifndef __MERGE_SETS__
#define __MERGE_SETS__

#include <vector>

class Merge_Sets{
public:
  // Create superset of n_element elements.
  // Also call a function call_on_merge when merging sets
  // with adams a1 and a2
  Merge_Sets(unsigned int n_element, void (*call_on_merge)
             (unsigned int, unsigned int));
  Merge_Sets(unsigned int n);

  // Destructor
  ~Merge_Sets();

  // Add more sets - (may be) useful for incremental algorithms
  void add_new_sets(unsigned n_new_element);
  
  // Print adams of all elements.
  void print_adams(void);
  
  // Merge sets containing elements e1 and e2.
  void merge(unsigned int e1, unsigned int e2);

  // Get adam for element e - performs update for faster queries in future
  unsigned int get_adam(unsigned int e);
  
  // Update adams for all or specific element e
  void update_adams();
  void update_adams(unsigned int e);
  
  // Private data
  unsigned int n_element;
  std::vector<unsigned int> adam;
  void (*call_on_merge)(unsigned int, unsigned int);
};

#endif
