// base class for merge criteria for agglomerative segmentation
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI.
//

#ifndef __IPUJ_AGGLO_MERGE_CRITERION__
#define __IPUJ_AGGLO_MERGE_CRITERION__

#include <ext/hash_map>

#include <segment.h>
#include <segment_agglo.h>

namespace std{
  using namespace __gnu_cxx;
}

namespace iput
{
namespace seg
{
  typedef std::vector<std::pair<Label,Label> > Adjacency_List;
  
  class Agglo_Merge_Criterion {
  public:
    virtual ~Agglo_Merge_Criterion(){}
    
    // functions to be declared by derived classes
    virtual bool 
      initialize_adjacency_statistics(const Array<unsigned char> * boundary_map,
				      const Array<Label> * segment_map) = 0;
    virtual double get_merge_priority(Label v_i, Label v_j) = 0;
    virtual bool should_be_merged(Label v_i, Label v_j, 
				  Merge_Priority f_threshold) = 0;
    virtual void unite_pairwise_statistics(Label n, Label source,
                                         Label target) = 0;
    virtual void move_pairwise_statistics(Label n, Label source,
                                         Label target) = 0;
    virtual void unite_unary_statistics(Label source, Label target) = 0;
    
    Adjacency_List adj_list;
    int is_verbose;
  };
}
}
#endif //__IPU_AGGLO_MERGE_CRITERION__

