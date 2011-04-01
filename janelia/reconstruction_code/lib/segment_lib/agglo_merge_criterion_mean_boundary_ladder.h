// merge criteria for agglomerative segmentation: mean boundary value
// with ladder (minimum area constraint)
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI.
//

#ifndef __IPUJ_AGGLO_MERGE_CRITERION_MEAN_BOUNDARY_LADDER__
#define __IPUJ_AGGLO_MERGE_CRITERION_MEAN_BOUNDARY_LADDER__

#include <ext/hash_map>

#include <agglo_merge_criterion.h>

namespace std{
  using namespace __gnu_cxx;
}

namespace iput
{
namespace seg
{
  class Agglo_Merge_Criterion_Mean_Boundary_Ladder:
    public Agglo_Merge_Criterion {
  public:
    Agglo_Merge_Criterion_Mean_Boundary_Ladder(){}
    ~Agglo_Merge_Criterion_Mean_Boundary_Ladder(){}
    
    // functions of base class being declared here
    bool initialize_adjacency_statistics(
      const Array<Boundary_Value> * boundary_map,
      const Array<Label> * segment_map);
    Merge_Priority get_merge_priority(Label v_i, Label v_j);
    bool should_be_merged(Label v_i, Label v_j, Merge_Priority f_threshold);
    void unite_pairwise_statistics(Label n, Label source, Label target);
    void move_pairwise_statistics(Label n, Label source, Label target);
    void unite_unary_statistics(Label source, Label target);

    unsigned int area_threshold;
  private:
    bool initialize_adjacency_statistics_2D__(
      const Array<Boundary_Value> * boundary_map,
      const Array<Label> * segment_map);
    bool initialize_adjacency_statistics_3D__(
      const Array<Boundary_Value> * boundary_map,
      const Array<Label> * segment_map);

    // boundary statistics
    Hash_Label_Pair_UInt32 segment_areas, sum_boundary_values,
      n_boundary_values;
  };
}
}
#endif //__IPU_AGGLO_MERGE_CRITERION_MEAN_BOUNDARY_LADDER__
