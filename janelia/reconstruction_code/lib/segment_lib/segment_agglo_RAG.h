// base class for agglomerative segmentation based on region adjacency graph
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI.
//

#ifndef __IPUJ_SEGMENT_AGGLO_RAG__
#define __IPUJ_SEGMENT_AGGLO_RAG__

#include <ext/hash_map>
#include <utility>                       // for std::pair
#include <algorithm>                     // for std::for_each
#include <boost/utility.hpp>             // for boost::tie
#include <boost/graph/graph_traits.hpp>  // for boost::graph_traits
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#include <hash_functions.h>

#include <segment.h>
#include <segment_agglo.h>
#include <agglo_merge_criterion.h>

namespace std{
  using namespace __gnu_cxx;
}

namespace iput
{
namespace seg
{
  // type of the graph: unsorted sequence of edges, unsorted set of vertices,
  // undirected graph.
  struct RAG_Vertex_Property{
    Label label;
  };
  typedef boost::adjacency_list<boost::vecS, boost::vecS,
                                boost::undirectedS, RAG_Vertex_Property> RAG;
  
  class Segment_Agglo_RAG: public Segment {
  public:
    Segment_Agglo_RAG(){
      segment_mappings = NULL;
    }
    ~Segment_Agglo_RAG(){
      if(segment_mappings != NULL)
        delete [] segment_mappings;
    }
    
    // functions to be declared by derived classes
    Agglo_Merge_Criterion * merge_criterion;

    // optional functions for profiling
    virtual void profile_begin_f_threshold(void){
#ifdef __DEBUG__
      std::cout << "Begin f threshold\n";
#endif
    }
    virtual void profile_end_f_threshold(void){
#ifdef __DEBUG__
      std::cout << "End f threshold\n";
#endif
    }
    virtual void profile_extracted_top_of_queue(void){
#ifdef __DEBUG__
      std::cout << "extracted top of queue\n";
#endif
    }
    virtual void profile_begin_merger(void){
#ifdef __DEBUG__
      std::cout << "Begin merger\n";
#endif
    }
    virtual void profile_end_merger(void){
#ifdef __DEBUG__
      std::cout << "End merger\n";
#endif
    }

    // functions of parent class that are being declared here
    bool initialize(const Array<Boundary_Value> * boundary_map, 
		    const Array<Label> * initial_segment_map);
    bool compute_segmentation(void);
    bool compute_segmentation_seeded(void);

    // variables to set by caller
    //
    // sequence of thresholds for which segmentation must be
    // saved. Must be in ascending order
    std::vector<Merge_Priority> f_thresholds;

    // Whether the mergers should cease after last f_threshold is
    // encountered
    bool is_enabled_break_at_last_f_threshold;
    
    // In case of seeded segmentation, list of segments that are
    // seeds.
    std::vector<Label> seeds;

    // variables having the results
    //
    // mapping of input segment labels to final segmentation. Number
    // of segment mappings is same as number of f_thresholds.
    Array<unsigned int> * segment_mappings;
    
  private:
    // merge priority queue
    Label_Pair_Q merge_q;
    std::pair<Label_Pair_Q_It, bool> insert_into_merge_q(
      Label v_i, Label v_j, Merge_Priority f);

    // the Region Adjacency Graph (RAG)
    RAG rag;
    typedef std::hash_map<Label, boost::graph_traits<RAG>::vertex_descriptor, 
      std::hash<Label>, eqlp> Hash_Label_Vertex;
    Hash_Label_Vertex label_to_vertex;
    

    // hash table for maintaining label-regions' adjacencies.
    Hash_Label_Pair_UInt32 label_adjacency;

    bool add_edge_into_RAG(Label v_i, Label v_j);

  };
}
}
#endif //__IPU_SEGMENT_AGGLO_RAG__
