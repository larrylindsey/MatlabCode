// Given a linkage graph generate segment-to-body-map. Negative
// segment ids are for dummy segments not to be included in the
// segment-to-body-map.
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI.
//

#include <mex.h>
#include <math.h>
#include <algorithm>
#include <utility>
#include <vector>
#include <ext/hash_map>
#include <map>
#include <merge_sets_h.h>
#include <boost/utility.hpp>             // for boost::tie
#include <boost/graph/graph_traits.hpp>  // for boost::graph_traits
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

using namespace __gnu_cxx;

struct equi{
  bool operator()(const unsigned int s1, const unsigned int s2) const
    {
      return s1==s2;
    }
};
typedef hash_map<unsigned int, unsigned int,
                 hash<unsigned int>,
                      equi> Hash_UInt32_UInt32;

struct equli{
  bool operator()(const unsigned long int s1, const unsigned long int s2) const
    {
      return s1==s2;
    }
};
typedef hash_map<unsigned long int, unsigned int,
                 hash<unsigned long int>,
                      equi> Hash_UInt64_UInt32;

typedef unsigned long int Label_Pair;
Label_Pair label_pair_2_id(unsigned int s1, unsigned int s2){
  if(s1>s2){
    Label_Pair lp1;
    lp1 = s2;
    lp1 <<= 32;
    Label_Pair lp;
    lp = s1;
    lp |= lp1;
    return lp;
  }
  Label_Pair lp1;
  lp1 = s1;
  lp1 <<= 32;
  Label_Pair lp;
  lp = s2;
  lp |= lp1;
  return lp;
}
std::pair<unsigned int, unsigned int> id_2_label_pair(Label_Pair lp){
  std::pair<unsigned int, unsigned int> s;
  s.first = lp >> 32;
  s.second = lp & 0x00000000ffffffff;
  return s;
}

typedef pair<unsigned long int, unsigned long int> Sec_Label_Pair;

// Priority queue for merging
struct ltdb{ bool operator()(const double f1, const double f2){ return f1<f2; }};
typedef std::map<const double, Label_Pair, ltdb> Label_Pair_Q;

// type of the graph: unsorted sequence of edges, unsorted set of vertices,
// undirected graph, no vertex property and edge property Boundary_Stat
using namespace boost;
typedef struct{
  double linkage_conf;
} Boundary_Stat_0;
enum edge_boundary_t {edge_boundary};
namespace boost{
  BOOST_INSTALL_PROPERTY(edge, boundary);
}
typedef property<edge_boundary_t, Boundary_Stat_0> Boundary_Stat;
typedef adjacency_list<setS, vecS, undirectedS, no_property, Boundary_Stat> Graph;

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    mexPrintf("input params:\n");
    mexPrintf("\t1.	linkage graph\n");
    mexPrintf("\t2. linkage threshold\n");
    mexPrintf("\t3. seed segments\n");
    mexPrintf("\t4. linkage threshold for seeding\n");
    mexPrintf("\tOutput:");
    mexPrintf("\t1. segment-to-body map Rx1 uint32\n");
    return;
  }
  if(nrhs!=4)
  {
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

  if(nlhs>1)
  {
    mexErrMsgTxt("Wrong number of outputs\n");
    return;
  }

  mexPrintf("START: get_sec_seg_2_body_map_from_links_3D_with_dummy_seeded_seg\n");
  const mxArray * linkage_graph_mx = prhs[0];
  const mxArray * linkage_threshold_mx = prhs[1];
  const mxArray * seeded_segs_mx = prhs[2];
  const mxArray * linkage_threshold_seed_mx = prhs[3];

  const int * size_1 = mxGetDimensions(linkage_graph_mx);
  int n_plane = std::max(size_1[0], size_1[1]);
  double linkage_threshold = * mxGetPr(linkage_threshold_mx);
  double linkage_threshold_seed = * mxGetPr(linkage_threshold_seed_mx);
  mexPrintf("linkage_threshold: %g, linkage_threshold_seed: %g\n", 
	    linkage_threshold, linkage_threshold_seed);

  const int section_bit_offset = 32;
  const unsigned long int segment_mask = 0x00000000ffffffff;

  mexPrintf("Build Merge_Sets_H datastructure for [section_id, segment_id] and\n");
  mexPrintf("and construct RAG\n");
  Merge_Sets_H<unsigned long int, hash<unsigned long int>, equli> M(NULL);
  {
    double * section_pair_links;
    int i;
    unsigned long int l1, l2, l;
    for(i=0; i<n_plane; i++){
      mexPrintf("plane: %d\n", i);
      section_pair_links = mxGetPr(mxGetCell(linkage_graph_mx, i));
      const int * size_l = mxGetDimensions(mxGetCell(linkage_graph_mx, i));
      int n_pair = size_l[0];
      mexPrintf("n_pair: %d\n", n_pair);
      int j;
      for(j=0; j<n_pair; j++){
        l1 = i+1;
        l1 <<= section_bit_offset;
        l = (unsigned int) (int) section_pair_links[j];
        l1 |= l;
        M.add_new_set_inc(l1);

        l2 = i+2;
        l2 <<= section_bit_offset;
        l = (unsigned int) (int) section_pair_links[j+n_pair];
        l2 |= l;
        M.add_new_set_inc(l2);

        if(section_pair_links[j]!=0 &&
           section_pair_links[j+n_pair]!=0 &&
           section_pair_links[j+2*n_pair]>linkage_threshold){
          M.merge(l1, l2);
        }
        
        // if dummy node then add an extra node in the next section in
        // case the body disappears for two sections.
        if(section_pair_links[j]!=0 &&
           section_pair_links[j+n_pair]<0){
	  l2 = i+3;
	  l2 <<= section_bit_offset;
	  l = (unsigned int) (int) section_pair_links[j+n_pair];
	  l2 |= l;
	  M.add_new_set_inc(l2);
	    
	  if(section_pair_links[j+2*n_pair]>linkage_threshold){
	    M.merge(l1, l2);
	  }
        }
      }
    }
  }

#ifdef __DEBUG__
  {
    mexPrintf("Body ids after linkage threshold:\n");
    Merge_Sets_H<unsigned long int, hash<unsigned long int>, equli>::
      Adam_Hash_Iterator it;
    unsigned long int l, l1;
    unsigned int s, p;
    for(it=M.adam.begin(); it!=M.adam.end(); it++){
      l = (*it).first;
      l1 = l & segment_mask;
      s = (int) l1;
      if(s<=0)
        continue;
      l1 = l >> section_bit_offset;
      p = (unsigned int) l1;
      mexPrintf("section: %u, segment: %u, body: %u\n", p, s, M.get_adam_id(l));
    }
  }
#endif

  mexPrintf("Marking bodies containing seeded segments.\n");
  Hash_UInt32_UInt32 is_seeded_body;
  {
    int i, j;
    unsigned long int l1, l;
    unsigned int * seeded_segs;
    for(i=0; i<n_plane; i++){
      seeded_segs = (unsigned int *) mxGetPr(mxGetCell(seeded_segs_mx, i));
      const int * size_l = mxGetDimensions(mxGetCell(seeded_segs_mx, i));
      int n_seed = std::max(size_l[0],size_l[1]);
      for(j=0; j<n_seed; j++){
	l1 = (unsigned long int) i+1;
	l1 <<= section_bit_offset;
	l = (unsigned long int) seeded_segs[j];
	l1 |= l;
	is_seeded_body[M.get_adam_id(l1)] = 1;
#ifdef __DEBUG__
	mexPrintf("is_seeded_body[%d] = 1\n", M.get_adam_id(l1));
#endif
      }
    }
  }

  mexPrintf("Constructing an undirected weighted graph of bodies with\n");
  mexPrintf("edge weights as linkage confidences. Adding all neighbors of\n");
  mexPrintf("seeded bodies to the merge q\n");
  Graph g(10000);
  property_map<Graph, edge_boundary_t>::type boundary_stats = get(edge_boundary_t(), g);
  Label_Pair_Q merge_q;
  Merge_Sets_H<unsigned int, hash<unsigned int>, equi> N(NULL);
  {
    double * section_pair_links;
    int i;
    unsigned long int l1, l2, l;
    for(i=0; i<n_plane; i++){
      mexPrintf("plane: %d\n", i);
      section_pair_links = mxGetPr(mxGetCell(linkage_graph_mx, i));
      const int * size_l = mxGetDimensions(mxGetCell(linkage_graph_mx, i));
      int n_pair = size_l[0];
      mexPrintf("n_pair: %d\n", n_pair);
      int j;
      for(j=0; j<n_pair; j++){
        l1 = i+1;
        l1 <<= section_bit_offset;
        l = (unsigned int) (int) section_pair_links[j];
        l1 |= l;

        l2 = i+2;
        l2 <<= section_bit_offset;
        l = (unsigned int) (int) section_pair_links[j+n_pair];
        l2 |= l;
	
	{
	  unsigned int s1 = M.get_adam_id(l1);
	  unsigned int s2 = M.get_adam_id(l2);
	  N.add_new_set_inc(s1);
	  N.add_new_set_inc(s2);
	  graph_traits<Graph>::edge_descriptor e;
	  bool is_inserted;
	  tie(e, is_inserted) = add_edge(s1, s2, g);
	  boundary_stats[e].linkage_conf = std::max(boundary_stats[e].linkage_conf,
					       section_pair_links[j+2*n_pair]);
#ifdef __DEBUG__
	  mexPrintf("edge <%u, %u, %g>\n", s1, s2, boundary_stats[e].linkage_conf);
#endif
	  if(is_seeded_body.find(s1)!=is_seeded_body.end() ||
	     is_seeded_body.find(s2)!=is_seeded_body.end()){
#ifdef __DEBUG__
	    mexPrintf("added to merge q\n");
#endif
	    // insert in merge queue
	    double m = -section_pair_links[j+2*n_pair];
	    Label_Pair_Q::iterator q_it = merge_q.find(m);
	    // keep priority keys unique.
	    while(q_it!=merge_q.end()){
	      m += 0.00001*((double)(rand()%1000)); // epsilon increments
	      q_it = merge_q.find(m);
	    }
	    {
	      Label_Pair lp = label_pair_2_id(s1, s2);
	      const std::pair<double, Label_Pair> v(m, lp);
	      merge_q.insert(v);
	    }
	  }
	}
      }
    }
  }

  mexPrintf("Merging bodies that do not include a seeded segment with those includeding seeded segments\n");
  {
    Label_Pair_Q::iterator q_it, q_it1;
    double m, m1; // merge criterion for a pair
    unsigned int s0, s1; // segments to be merged
    unsigned int n0, n1; // their respective neighbors
    graph_traits<Graph>::out_edge_iterator edge_0, edge_0_end,
      edge_1, edge_1_end; // to iterate over the neighbors of s0, s1
    graph_traits<Graph>::edge_descriptor e0, e1;
    unsigned int b, n; // boundary value and number of pixels

    // get number of neighbors for each vertex
    Hash_UInt32_UInt32 n_neighbor;
    {
      graph_traits<Graph>::vertex_iterator vs, ve;
      tie(vs, ve) = vertices(g);
      for(; vs!=ve; vs++)
        n_neighbor[(*vs)] = out_degree(*vs, g);
    }
    q_it = merge_q.begin();
    if(q_it!=merge_q.end()){
      m = -(*q_it).first;
      tie(s0,s1) = id_2_label_pair((*q_it).second);
#ifdef __DEBUG__
      mexPrintf("m: %g, s0: %u, s1: %u\n", m, s0, s1);
#endif
      while(q_it!=merge_q.end() && m>linkage_threshold_seed){
	// merge segments s0 and s1
	if(s0==0 || s1==0 || (is_seeded_body[s0]>0 && is_seeded_body[s1]>0) ||
	   (N.get_adam_id(s0)==N.get_adam_id(s1))){
#ifdef __DEBUG__
	  mexPrintf("skipping.\n");
#endif
	  merge_q.erase(q_it);
	  q_it = merge_q.begin();
	  if(q_it==merge_q.end())
	    break;
	  m = -(*q_it).first;
	  tie(s0,s1) = id_2_label_pair((*q_it).second);
#ifdef __DEBUG__
	  mexPrintf("m: %g, s0: %u, s1: %u\n", m, s0, s1);
#endif
	  continue;
	}
        
#ifdef __DEBUG__
	mexPrintf("merging\n");
#endif
	// Exactly one of s0 and s1 is seeded. If s1 is seeded then swap them.
	if(is_seeded_body[s1]>0){
	  unsigned int s01 = s0;
	  s0 = s1;
	  s1 = s01;
	}
	// From now on, s0 is seeded and s1 is not seeded.
        
	// merge the body of s1 into s0
	N.merge(s0,s1);

	// remove edge from graph
	remove_edge(s0, s1, g);
        
	// remove from merge queue
	merge_q.erase(q_it);

	// set the merged body to be seeded.
	is_seeded_body[s1]=1;

	// add the neighbors of s1 to the merge q.
	tie(edge_1, edge_1_end) = out_edges(s1, g);
	while(edge_1!=edge_1_end){
	  n1 = target(*edge_1, g);

	  double m = -boundary_stats[*edge_1].linkage_conf;
	  Label_Pair_Q::iterator q_it = merge_q.find(m);
	  // keep priority keys unique.
	  while(q_it!=merge_q.end()){
	    m += 0.00001*((double)(rand()%1000)); // epsilon increments
	    q_it = merge_q.find(m);
	  }
	  {
	    Label_Pair lp = label_pair_2_id(s1, n1);
	    const std::pair<double, Label_Pair> v(m, lp);
	    merge_q.insert(v);
	  }
	  edge_1++;
	}

	q_it = merge_q.begin();
	if(q_it==merge_q.end())
	  break;
	m = -(*q_it).first;
	tie(s0,s1) = id_2_label_pair((*q_it).second);
#ifdef __DEBUG__
	  mexPrintf("m: %g, s0: %u, s1: %u\n", m, s0, s1);
#endif
      }
    }
  }

  mexPrintf("Construct segment-to-body-map\n");
  {
    mexPrintf("get the number of valid segments and bodies.\n");
    unsigned int n_segment_real = 0;
    Hash_UInt32_UInt32 body_id_hash;
    unsigned int a, n, p;
    unsigned long int l, l1;
    int s;
    Merge_Sets_H<unsigned long int, hash<unsigned long int>, equli>::
      Adam_Hash_Iterator it;
    for(it=M.adam.begin(); it!=M.adam.end(); it++){
      l = (*it).first;
      l1 = l & segment_mask;
      s = (int) l1;
      if(s<=0)
        continue;
      l1 = l >> section_bit_offset;
      p = (unsigned int) l1;
#ifdef __DEBUG__
      mexPrintf("section: %u, segment: %u\n", p, s);
#endif
      n_segment_real++;
      a = N.get_adam_id(M.get_adam_id(l));
      if(body_id_hash[a]==0){
        n = body_id_hash.size();
        body_id_hash[a] = n;
      }
    }

    plhs[0] = mxCreateNumericMatrix(n_segment_real, 3, mxUINT32_CLASS, mxREAL);
    unsigned int * sec_seg_2_body = (unsigned int *) mxGetPr(plhs[0]);
    for(it=M.adam.begin(), n=0; it!=M.adam.end(); it++){
      l = (*it).first;
      l1 = l & segment_mask;
      s = (int) l1;
      if(s<=0)
        continue;
      l1 = l >> section_bit_offset;
      p = (unsigned int) l1;
      sec_seg_2_body[n] = p; // plane or layer id
      sec_seg_2_body[n+n_segment_real] = s; // segment id
      sec_seg_2_body[n+2*n_segment_real] = body_id_hash[N.get_adam_id(M.get_adam_id(l))];
      
      n++;
    }
  }
  mexPrintf("STOP: get_sec_seg_2_body_map_from_links_3D_with_dummy_seeded_seg\n");
  return;
}
