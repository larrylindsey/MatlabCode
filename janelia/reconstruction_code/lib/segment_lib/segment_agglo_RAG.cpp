// base class for agglomerative segmentation based on region adjacency graph
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI.
//

#include <utility>                       // for std::pair
#include <algorithm>                     // for std::for_each
#include <boost/utility.hpp>             // for boost::tie
#include <boost/graph/graph_traits.hpp>  // for boost::graph_traits
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#include <hash_functions.h>
#include <merge_sets_h.h>

#include <segment_agglo_RAG.h>
#include <agglo_merge_criterion.h>

namespace iput
{
  namespace seg
  {

    bool Segment_Agglo_RAG::initialize(
      const Array<Boundary_Value> * boundary_map,
      const Array<Label> * initial_segment_map){
      bool is_valid = merge_criterion->
        initialize_adjacency_statistics(boundary_map, initial_segment_map);
      if(is_valid){
        Adjacency_List::iterator vp;
        for(vp = merge_criterion->adj_list.begin(); 
            vp != merge_criterion->adj_list.end(); vp++)
          if(!add_edge_into_RAG((*vp).first, (*vp).second))
            return false;
      }
      return is_valid;
    }

    bool Segment_Agglo_RAG::add_edge_into_RAG(Label v_i, Label v_j){
      if(label_to_vertex.find(v_i)==label_to_vertex.end()){
        RAG_Vertex_Property p;
        p.label = v_i;
        label_to_vertex[v_i] = boost::add_vertex(p, rag);
      }
      if(label_to_vertex.find(v_j)==label_to_vertex.end()){
        RAG_Vertex_Property p;
        p.label = v_j;
        label_to_vertex[v_j] = boost::add_vertex(p, rag);
      }

      if(v_i==v_j || v_i==0 || v_j==0)
        return true;
      
      boost::graph_traits<RAG>::edge_descriptor e;
      bool is_inserted;
      boost::tie(e, is_inserted) = boost::add_edge(label_to_vertex[v_i], 
                                                   label_to_vertex[v_j], rag);
      return is_inserted;
    }
    
  
    bool Segment_Agglo_RAG::compute_segmentation(void){

      segment_mappings = (Array<unsigned int> *) new
        Array<unsigned int>[f_thresholds.size()+1];
      
      Merge_Sets_H<unsigned int, std::hash<unsigned int>, equi>
        segment_sets(NULL);
      Hash_Label_Pair_Q_It merge_q_pos;
      // build merge queue
      {
        boost::graph_traits<RAG>::edge_iterator e_s, e_t;
        boost::tie(e_s, e_t) = edges(rag);
        for(; e_s!=e_t; e_s++){
          Label v_i, v_j;
          v_i = rag[boost::source(*e_s, rag)].label;
          v_j = rag[boost::target(*e_s, rag)].label;
          Merge_Priority m = merge_criterion->get_merge_priority(v_i, v_j);
          Label_Pair_Q_It q_it;
          bool is_inserted;
          boost::tie(q_it,is_inserted) = insert_into_merge_q(v_i, v_j, m);
          if(!is_inserted)
            return false;
          segment_sets.add_new_set_inc(v_i);
          segment_sets.add_new_set_inc(v_j);
          Label_Pair_Key lpk = lp_2_lpk(Label_Pair(v_i, v_j));
          label_adjacency[lpk] = 1;
          merge_q_pos[lpk] = q_it;
        }
      }

      // make mergers. Save label mapping for f_threshold_seq.
      {
        Label_Pair_Q::iterator q_it, q_it1;
        Merge_Priority m, m1; // merge criterion for a pair
        Label s0, s1; // segments to be merged
        Label n0, n1; // their respective neighbors
        int f_threshold_id=0;
        double f_threshold=0;
        boost::graph_traits<RAG>::out_edge_iterator edge_1, edge_1_end;
        // to iterate over the neighbors of s0, s1
        boost::graph_traits<RAG>::edge_descriptor e0, e1;

        // get number of neighbors for each vertex
        Hash_Label_UInt32 n_neighbor;
        {
          boost::graph_traits<RAG>::vertex_iterator vs, ve;
          boost::tie(vs, ve) = vertices(rag);
          for(; vs!=ve; vs++)
            n_neighbor[rag[(*vs)].label] = boost::out_degree(*vs, rag);
        }
        q_it = merge_q.begin();
        while(q_it!=merge_q.end()){
          if(f_threshold_id<f_thresholds.size())
            f_threshold = f_thresholds[f_threshold_id];
	  if(is_verbose>16)
	    std::cout << "f_threshold: " << f_threshold << '\n';
          m = (*q_it).first;
          boost::tie(s0,s1) = (*q_it).second;
          merge_q.erase(q_it);
          q_it = merge_q.begin();

          if(s0==0 || s1==0 || s0==s1 || 
             merge_criterion->should_be_merged(s0, s1, f_threshold)){
            profile_begin_merger();
            
            if(s0==0 || s1==0 || s0==s1 ||
               label_adjacency[lp_2_lpk(Label_Pair(s0,s1))]==0)
              continue;
        
            // merge the superpixel sets of s1 into s0
            segment_sets.merge(s0,s1);

            // merge the neighbors of segments s0 and s1
            if(n_neighbor[s0]<n_neighbor[s1]){
              Label s01 = s0;
              s0 = s1;
              s1 = s01;
            }
            n_neighbor[s0] += n_neighbor[s1];

            merge_criterion->unite_unary_statistics(s1, s0);
            
	    if(is_verbose>16)
	      std::cout << "s0: " << s0 << " s1: " << s1 << '\n';
            tie(edge_1, edge_1_end) = out_edges(label_to_vertex[s1], rag);
            while(edge_1!=edge_1_end){
              n1 = rag[target(*edge_1, rag)].label;
	      if(is_verbose>16)
		std::cout << "n1: " << n1 << '\n';
              Label_Pair_Key lpk_n1_s1 = lp_2_lpk(Label_Pair(n1,s1));
              if(n1==s0 || label_adjacency[lpk_n1_s1]==0){
                edge_1++;
                continue;
              }
            
              if(n_neighbor[n1]!=0){
		if(is_verbose>16)
		  std::cout << "n1 still valid\n";

                // n1 is valid and is a neighbor of s1
                Label_Pair_Key lpk_n1_s0 = lp_2_lpk(Label_Pair(n1,s0));
                if(label_adjacency[lpk_n1_s0]==1){
		  if(is_verbose>16)
		    std::cout << "merging stats of " << n1 << "-" << s1
			      << " into " << n1 << "-" << s0 << '\n';

                  // n1 is a neighbor of s0. Remove adjacency between n1
                  // and s1 and merge the statistics of n1 with s0 and
                  // s1 into s0

                  // invalidate the merge queue positions of n1-s1 and
                  // n1-s0
                  (*merge_q_pos[lpk_n1_s1]).second.first = 0;
                  (*merge_q_pos[lpk_n1_s1]).second.second = 0;
                  (*merge_q_pos[lpk_n1_s0]).second.first = 0;
                  (*merge_q_pos[lpk_n1_s0]).second.second = 0;
                  
                  label_adjacency[lpk_n1_s1] = 0;

                  merge_criterion->unite_pairwise_statistics(n1, s1, s0);
                
                  // insert in merge queue
                  m1 = merge_criterion->get_merge_priority(n1, s0);
                  Label_Pair_Q_It q_it;
                  bool is_inserted;
                  boost::tie(q_it,is_inserted) =
                    insert_into_merge_q(n1, s0, m1);
                  if(!is_inserted)
                    return false;
                  merge_q_pos[lpk_n1_s0] = q_it;
                }
                else{
                  
		  if(is_verbose>16)
		    std::cout << "moving stats of " << n1 << '-' << s1
			      << " to " << n1 << '-' << s0 <<'\n';

                  // move this neighbor of s1 to s0
                  boost::graph_traits<RAG>::edge_descriptor e;
                  bool is_inserted;
                  boost::tie(e, is_inserted) = boost::add_edge(
                    label_to_vertex[s0], label_to_vertex[n1], rag);

                  label_adjacency[lpk_n1_s1] = 0;
                  
                  label_adjacency[lpk_n1_s0] = 1;
                  merge_criterion->move_pairwise_statistics(n1, s1, s0);

                  merge_q_pos[lpk_n1_s0] = merge_q_pos[lpk_n1_s1];
                  (*merge_q_pos[lpk_n1_s0]).second.first = n1;
                  (*merge_q_pos[lpk_n1_s0]).second.second = s0;
                }
              }
          
              edge_1++;
            }

            // invalidate s1
            n_neighbor[s1] = 0;
        
            profile_end_merger();
          }

          // remove edge from graph
          label_adjacency[lp_2_lpk(Label_Pair(s0,s1))]=0;
        

          if(m>=f_threshold){
            if(f_threshold_id>=f_thresholds.size() &&
               is_enabled_break_at_last_f_threshold)
              break;

            if(f_threshold_id<f_thresholds.size()){
              // copy the current superpixel-to-segment mapping into
              // output variable
	      if(is_verbose>2)
		std::cout << "Copy for f_threshold: " << f_threshold << '\n';
              segment_mappings[f_threshold_id].n_dimension = 2;
              segment_mappings[f_threshold_id].dimensions.push_back(2);
              segment_mappings[f_threshold_id].dimensions.push_back
		(n_neighbor.size());
              segment_mappings[f_threshold_id].buffer = (unsigned int*)
                new unsigned int[2*n_neighbor.size()];
              Hash_Label_UInt32::iterator hi;
              int i;
              for(i=0, hi=n_neighbor.begin(); hi!=n_neighbor.end();
                  hi++, i++){
                segment_mappings[f_threshold_id].buffer[2*i] = (*hi).first;
                segment_mappings[f_threshold_id].buffer[2*i+1] =
                  segment_sets.get_adam((*hi).first);
              }
            }
            
            f_threshold_id++;
          }
        }

        if(!is_enabled_break_at_last_f_threshold && q_it==merge_q.end()){
          // copy last superpixel-to-segment mapping into output variable
	  if(is_verbose>2)
	    std::cout << "Copy for end segmentation\n";
          int j = f_thresholds.size();
          segment_mappings[j].n_dimension = 2;
          segment_mappings[j].dimensions.push_back(2);
          segment_mappings[j].dimensions.push_back(n_neighbor.size());
          segment_mappings[j].buffer = (unsigned int*)
            new unsigned int[2*n_neighbor.size()];
          Hash_Label_UInt32::iterator hi;
          int i;
          for(i=0, hi=n_neighbor.begin(); hi!=n_neighbor.end();
              hi++, i++){
            segment_mappings[j].buffer[2*i] = (*hi).first;
            segment_mappings[j].buffer[2*i+1] =
              segment_sets.get_adam((*hi).first);
          }
        }
      }
      return true;
    }

    bool Segment_Agglo_RAG::compute_segmentation_seeded(void){

      if(seeds.empty())
	return false;
      // create hash for seeds
      Hash_Label_UInt32 is_seeded;
      {
	std::vector<Label>::iterator i;
	for(i=seeds.begin(); i!=seeds.end(); i++)
	  is_seeded[*i]=1;
      }

      segment_mappings = (Array<unsigned int> *) new
        Array<unsigned int>[f_thresholds.size()+1];
      
      Merge_Sets_H<unsigned int, std::hash<unsigned int>, equi>
        segment_sets(NULL);
      Hash_Label_Pair_Q_It merge_q_pos;
      Hash_Label_UInt32 n_neighbor;
      // initialize merge_set data-structure, label adjacency and
      // build merge queue
      {
        boost::graph_traits<RAG>::edge_iterator e_s, e_t;
        boost::tie(e_s, e_t) = edges(rag);
        for(; e_s!=e_t; e_s++){
          Label v_i, v_j;
          v_i = rag[boost::source(*e_s, rag)].label;
          v_j = rag[boost::target(*e_s, rag)].label;
	  segment_sets.add_new_set_inc(v_i);
	  segment_sets.add_new_set_inc(v_j);
	  Label_Pair_Key lpk = lp_2_lpk(Label_Pair(v_i, v_j));
	  label_adjacency[lpk] = 1;
	  if((is_seeded[v_i]==1&&is_seeded[v_j]!=1) || 
	     (is_seeded[v_i]!=1&&is_seeded[v_j]==1)){
	    if(is_verbose>16)
	      std::cout << "merge_q.insert("<<v_i<<","<<v_j<<")\n";
	    // add this pair only if exactly one of them is seeded
	    Merge_Priority m = merge_criterion->get_merge_priority(v_i, v_j);
	    Label_Pair_Q_It q_it;
	    bool is_inserted;
	    boost::tie(q_it,is_inserted) = insert_into_merge_q(v_i, v_j, m);
	    if(!is_inserted)
	      return false;
	    merge_q_pos[lpk] = q_it;
	  }
        }
      }

      // make mergers. Save label mapping for f_threshold_seq.
      {
        Label_Pair_Q::iterator q_it, q_it1;
        Merge_Priority m, m1; // merge criterion for a pair
        Label s0, s1; // segments to be merged
        Label n0, n1; // their respective neighbors
        int f_threshold_id=0;
        double f_threshold=0;
        boost::graph_traits<RAG>::out_edge_iterator edge_1, edge_1_end;
        // to iterate over the neighbors of s0, s1
        boost::graph_traits<RAG>::edge_descriptor e0, e1;

        // get number of neighbors for each vertex
        {
          boost::graph_traits<RAG>::vertex_iterator vs, ve;
          boost::tie(vs, ve) = vertices(rag);
          for(; vs!=ve; vs++)
            n_neighbor[rag[(*vs)].label] = boost::out_degree(*vs, rag);
        }
        q_it = merge_q.begin();
        while(q_it!=merge_q.end()){
          if(f_threshold_id<f_thresholds.size())
            f_threshold = f_thresholds[f_threshold_id];
	  if(is_verbose>16)
	    std::cout << "f_threshold: " << f_threshold << '\n';
          m = (*q_it).first;
          boost::tie(s0,s1) = (*q_it).second;
          merge_q.erase(q_it);
          q_it = merge_q.begin();

          if(s0==0 || s1==0 || s0==s1 || 
             merge_criterion->should_be_merged(s0, s1, f_threshold)){
            profile_begin_merger();
            
            if(s0==0 || s1==0 || s0==s1 ||
               label_adjacency[lp_2_lpk(Label_Pair(s0,s1))]==0)
              continue;
        
	    if((is_seeded[s0]==1&&is_seeded[s1]==1) ||
	       (is_seeded[s0]!=1&&is_seeded[s1]!=1)){
	      std::cout << "is_seeded[s0]: " << is_seeded[s0] << 
		" , is_seeded[s1]: " << is_seeded[s1] << '\n';
	      return false;
	    }

            // merge the superpixel sets of s1 into s0
            segment_sets.merge(s0,s1);

	    // ensure that after this s0 is seeded and s1 is not
            if(is_seeded[s1]){
              Label s01 = s0;
              s0 = s1;
              s1 = s01;
            }
            n_neighbor[s0] += n_neighbor[s1];

            merge_criterion->unite_unary_statistics(s1, s0);
            
            // merge the neighbors of segments s1 into s0
	    if(is_verbose>16)
	      std::cout << "s0: " << s0 << " s1: " << s1 << '\n';
            tie(edge_1, edge_1_end) = out_edges(label_to_vertex[s1], rag);
            while(edge_1!=edge_1_end){
              n1 = rag[target(*edge_1, rag)].label;
	      if(is_verbose>16)
		std::cout << "n1: " << n1 << '\n';
              Label_Pair_Key lpk_n1_s1 = lp_2_lpk(Label_Pair(n1,s1));
              if(n1==s0 || label_adjacency[lpk_n1_s1]==0){
                edge_1++;
                continue;
              }
            
              if(n_neighbor[n1]!=0){
		if(is_verbose>16)
		  std::cout << "n1 still valid\n";

                // n1 is valid and is a neighbor of s1
                Label_Pair_Key lpk_n1_s0 = lp_2_lpk(Label_Pair(n1,s0));
                if(label_adjacency[lpk_n1_s0]==1){
		  if(is_verbose>16)
		    std::cout << "merging stats of " << n1 << "-" << s1
			      << " into " << n1 << "-" << s0 << '\n';

                  // n1 is a neighbor of s0. Remove adjacency between n1
                  // and s1 and merge the statistics of n1 with s0 and
                  // s1 into s0

		  if(is_seeded[n1]!=1){
		    if(is_verbose>16)
		      std::cout << "invalidate the merge queue positions " << 
			"of n1("<<n1<<")-s0("<<s0<<")\n";
		    (*merge_q_pos[lpk_n1_s0]).second.first = 0;
		    (*merge_q_pos[lpk_n1_s0]).second.second = 0;
                  }
		  if(is_seeded[n1]==1){
		    if(is_verbose>16)
		      std::cout << "invalidate the merge queue positions " << 
			"of n1("<<n1<<")-s1("<<s1<<")\n";
		    (*merge_q_pos[lpk_n1_s1]).second.first = 0;
		    (*merge_q_pos[lpk_n1_s1]).second.second = 0;
		  }

                  label_adjacency[lpk_n1_s1] = 0;

		  if(is_verbose>16)
		    std::cout << "unite statistics of n1-s1 and n1-s0\n";
                  merge_criterion->unite_pairwise_statistics(n1, s1, s0);

		  if(!is_seeded[n1]){
		    if(is_verbose>16)
		      std::cout << "n1 not seeded so inserting in merge Q\n";
		    m1 = merge_criterion->get_merge_priority(n1, s0);
		    Label_Pair_Q_It q_it;
		    bool is_inserted;
		    boost::tie(q_it,is_inserted) =
		      insert_into_merge_q(n1, s0, m1);
		    if(!is_inserted)
		      return false;
		    merge_q_pos[lpk_n1_s0] = q_it;
		  }
                }
                else{
                  
		  if(is_verbose>16)
		    std::cout << "moving stats of " << n1 << '-' << s1
			      << " to " << n1 << '-' << s0 <<'\n';
		  // move this neighbor of s1 to s0
		  boost::graph_traits<RAG>::edge_descriptor e;
		  bool is_inserted;
		  boost::tie(e, is_inserted) = 
		    boost::add_edge(label_to_vertex[s0], 
				    label_to_vertex[n1], rag);

		  label_adjacency[lpk_n1_s1] = 0;
                  
		  label_adjacency[lpk_n1_s0] = 1;
		  merge_criterion->move_pairwise_statistics(n1, s1, s0);

		  if(is_seeded[n1]==1){
		    if(is_verbose>16)
		      std::cout << "n1 seed. Removing n1-s1 from merge Q.\n";
		    (*merge_q_pos[lpk_n1_s1]).second.first = 0;
		    (*merge_q_pos[lpk_n1_s1]).second.second = 0;
		  }
		  if(is_seeded[n1]!=1){
		    if(is_verbose>16)
		      std::cout << "n1 not seeded. Adding n1-s0 to merge Q.\n";
		    m1 = merge_criterion->get_merge_priority(n1, s0);
		    Label_Pair_Q_It q_it;
		    bool is_inserted;
		    boost::tie(q_it,is_inserted) =
		      insert_into_merge_q(n1, s0, m1);
		    if(!is_inserted)
		      return false;
		    merge_q_pos[lpk_n1_s0] = q_it;
		  }
                }
              }
          
              edge_1++;
            }

            // invalidate s1
            n_neighbor[s1] = 0;
	    
            profile_end_merger();
          }

          // remove edge from graph
          label_adjacency[lp_2_lpk(Label_Pair(s0,s1))]=0;
        

          if(m>=f_threshold){
            if(f_threshold_id>=f_thresholds.size() &&
               is_enabled_break_at_last_f_threshold)
              break;

            if(f_threshold_id<f_thresholds.size()){
              // copy the current superpixel-to-segment mapping into
              // output variable
	      if(is_verbose>2)
		std::cout << "Copy for f_threshold: " << f_threshold << '\n';
              segment_mappings[f_threshold_id].n_dimension = 2;
              segment_mappings[f_threshold_id].dimensions.push_back(2);
              segment_mappings[f_threshold_id].dimensions.push_back
		(n_neighbor.size());
              segment_mappings[f_threshold_id].buffer = (unsigned int*)
                new unsigned int[2*n_neighbor.size()];
              Hash_Label_UInt32::iterator hi;
              int i;
              for(i=0, hi=n_neighbor.begin(); hi!=n_neighbor.end();
                  hi++, i++){
                segment_mappings[f_threshold_id].buffer[2*i] = (*hi).first;
                segment_mappings[f_threshold_id].buffer[2*i+1] =
                  segment_sets.get_adam((*hi).first);
              }
            }
            
            f_threshold_id++;
          }
        }

        if(!is_enabled_break_at_last_f_threshold && q_it==merge_q.end()){
          // copy last superpixel-to-segment mapping into output variable
	  if(is_verbose>2)
	    std::cout << "Copy for end segmentation\n";
          int j = f_thresholds.size();
          segment_mappings[j].n_dimension = 2;
          segment_mappings[j].dimensions.push_back(2);
          segment_mappings[j].dimensions.push_back(n_neighbor.size());
          segment_mappings[j].buffer = (unsigned int*)
            new unsigned int[2*n_neighbor.size()];
          Hash_Label_UInt32::iterator hi;
          int i;
          for(i=0, hi=n_neighbor.begin(); hi!=n_neighbor.end();
              hi++, i++){
            segment_mappings[j].buffer[2*i] = (*hi).first;
            segment_mappings[j].buffer[2*i+1] =
              segment_sets.get_adam((*hi).first);
          }
        }
      }
      return true;
    }
  
    std::pair<Label_Pair_Q_It, bool> Segment_Agglo_RAG::insert_into_merge_q(
      Label v_i, Label v_j, Merge_Priority f){
      Label_Pair_Q_It q_it = merge_q.find(f); // keep priority keys unique.
      while(q_it!=merge_q.end()){
        f += 0.00001*((double)(rand()%1000)); // epsilon increments
        q_it = merge_q.find(f);
      }

      bool is_inserted;
      {
        Label_Pair lp(v_i, v_j);
        const std::pair<double, Label_Pair> v(f, lp);
        boost::tie(q_it, is_inserted)  = merge_q.insert(v);
      }
      return std::pair<Label_Pair_Q_It,bool>(q_it, is_inserted);
    }
  }
}

