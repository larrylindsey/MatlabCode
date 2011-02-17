// merge criteria for agglomerative segmentation: mean boundary value
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI.
//

#include <iostream>

#include <segment.h>
#include <agglo_merge_criterion_mean_boundary_ladder.h>



#define PARENT_CLASS__ Agglo_Merge_Criterion_Mean_Boundary_Ladder

#define COLLECT_BOUNDARY_STATISTICS__ {                         \
    sum_boundary_values[LP_2_LPK(Label_Pair(s1,s2))] += *b_ptr; \
    n_boundary_values[LP_2_LPK(Label_Pair(s1,s2))] ++;          \
  }

#define COLLECT_INTERNAL_STATISTICS__ {         \
    segment_areas[s0]++;                        \
  }

#include <agglo_merge_criterion_initialize_adjacency_statistics_2D.h>
#include <agglo_merge_criterion_initialize_adjacency_statistics_3D.h>



namespace iput
{
  namespace seg
  {
    bool Agglo_Merge_Criterion_Mean_Boundary_Ladder::initialize_adjacency_statistics(
      const Array<Boundary_Value> * boundary_map,
      const Array<Label> * segment_map){

      if(boundary_map->n_dimension==2){
        bool is_valid = initialize_adjacency_statistics_2D__(boundary_map, 
                                                             segment_map);
        return is_valid;
      }
      if(boundary_map->n_dimension==3){
        bool is_valid = initialize_adjacency_statistics_3D__(boundary_map, 
                                                             segment_map);
        return is_valid;
      }

      return false;
    }
    
    Merge_Priority Agglo_Merge_Criterion_Mean_Boundary_Ladder::get_merge_priority(Label v_i, Label v_j){
      Label_Pair_Key lpk = LP_2_LPK(Label_Pair(v_i,v_j));
      unsigned int n = n_boundary_values[lpk];
      if(n!=0)
        return (Merge_Priority)((double)sum_boundary_values[lpk])/
	  ((double)n);
      else
        return (Merge_Priority) 100000;
    }
    
    bool Agglo_Merge_Criterion_Mean_Boundary_Ladder::should_be_merged(Label v_i, Label v_j, Merge_Priority f_threshold){
      return get_merge_priority(v_i, v_j)<f_threshold || std::min(segment_areas[v_i], segment_areas[v_j])<area_threshold;
    }
    
    void Agglo_Merge_Criterion_Mean_Boundary_Ladder::unite_pairwise_statistics(Label n, Label source, Label target){
      Label_Pair_Key lpk_ns = LP_2_LPK(Label_Pair(n,source));
      Label_Pair_Key lpk_nt = LP_2_LPK(Label_Pair(n,target));
      sum_boundary_values[lpk_nt] += sum_boundary_values[lpk_ns];
      n_boundary_values[lpk_nt] += n_boundary_values[lpk_ns];
    }
    
    void Agglo_Merge_Criterion_Mean_Boundary_Ladder::move_pairwise_statistics(Label n, Label source, Label target){
      Label_Pair_Key lpk_ns = LP_2_LPK(Label_Pair(n,source));
      Label_Pair_Key lpk_nt = LP_2_LPK(Label_Pair(n,target));
      sum_boundary_values[lpk_nt] = sum_boundary_values[lpk_ns];
      n_boundary_values[lpk_nt] = n_boundary_values[lpk_ns];
    }

    void Agglo_Merge_Criterion_Mean_Boundary_Ladder::unite_unary_statistics(Label source, Label target){
      segment_areas[target] += segment_areas[source];
#ifdef __DEBUG__
      std::cout << "segment_areas[" << target << "]: "
                << segment_areas[target] << '\n';
#endif
    }
  }  
}
