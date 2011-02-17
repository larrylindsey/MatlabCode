// base class for agglomerative segmentation based on region adjacency graph
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI.
//

#ifndef __IPUJ_SEGMENT_AGGLO__
#define __IPUJ_SEGMENT_AGGLO__

#include <ext/hash_map>
#include <utility>                       // for std::pair
#include <algorithm>                     // for std::for_each
#include <boost/utility.hpp>             // for boost::tie
#include <boost/graph/graph_traits.hpp>  // for boost::graph_traits
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#include <hash_functions.h>

namespace std{
  using namespace __gnu_cxx;
}

namespace iput
{
namespace seg
{
  typedef double Merge_Priority;
  struct ltmc{ bool operator()(const Merge_Priority f1,
                               const Merge_Priority f2){
    return f1<f2; }};
  
  // Priority queue of label pairs
  typedef std::map<const Merge_Priority, Label_Pair, ltmc> Label_Pair_Q;
  typedef Label_Pair_Q::iterator Label_Pair_Q_It;
  typedef std::hash_map<Label_Pair_Key, Label_Pair_Q_It, std::hash<Label_Pair_Key>, eqlpk> Hash_Label_Pair_Q_It;
}
}
#endif //__IPU_SEGMENT_AGGLO__
