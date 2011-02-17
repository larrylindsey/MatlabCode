#include <mex.h>
#include <math.h>
#include <stdio.h>
#include <set>
#include <ext/hash_map>
#include <utility>                       // for std::pair
#include <algorithm>                     // for std::for_each
#include <boost/utility.hpp>             // for boost::tie
#include <boost/graph/graph_traits.hpp>  // for boost::graph_traits
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <queue>

#include <merge_sets_h.h>

namespace std{
  using namespace __gnu_cxx;
}

using namespace boost;

unsigned long int get_triplet_word(unsigned int a, unsigned int b, unsigned int c){
  // ensure a <= b <= c
  if(a>b){
    unsigned int s = a;
    a = b;
    b = s;
  }
  if(b>c){
    unsigned int s = b;
    b = c;
    c = s;
  }
  if(a>b){
    unsigned int s = a;
    a = b;
    b = s;
  }
  
  unsigned long int l, m;
  l = (unsigned long int) a;
  l <<= 32;
  m = (unsigned long int) b;
  m <<= 16;
  l |= m;
  m = (unsigned long int) c;
  l |= m;
  return l;
}

typedef std::pair<unsigned int, unsigned int> Pair_UInt32;

//hash function 
struct equli{ bool operator()(const unsigned long int s1, const unsigned long int s2) const
    {return s1==s2;}};
typedef  std::hash_map<unsigned long int, unsigned int,
                       std::hash<unsigned long int>, equli> Hash_UInt64_UInt32;

// Hash function 
struct equi{ bool operator()(const unsigned int s1, const unsigned int s2) const
    {return s1==s2;}};
typedef  std::hash_map<unsigned int, unsigned int,
                       std::hash<unsigned int>, equi> Hash_UInt32_UInt32;
// Hash function 
typedef  std::hash_map<unsigned int, double,
                       std::hash<unsigned int>, equi> Hash_UInt32_Double;

// type of the graph: unsorted sequence of edges, unsorted set of vertices,
// undirected graph, no vertex property and edge property Boundary_Stat
typedef adjacency_list<vecS, vecS, undirectedS, no_property> Graph;

struct ltui{ bool operator()(const unsigned int s1, const unsigned int s2) const
    {return s1<s2;}};

/*
unsigned int get_adj_id(unsigned int i, unsigned int j,
                           unsigned int n_node){
  if(i<j)
    return i + j*n_node;
  else
    return j + i*n_node;
}
*/

unsigned int get_adj_id_32(unsigned int i, unsigned int j){
  if(i<j)
    return (i<<16) | j;
  else
    return (j<<16) | i;
}

class Weighted_Edge{
	public:
  double w;
  unsigned int i, j;
};
bool operator<(const Weighted_Edge& a, const Weighted_Edge& b){
	return a.w<b.w;
}


void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    mexPrintf("input params:\n");
    mexPrintf("\t1.	NxN adjacency matrix\n");
    return;
  }
  if(nrhs!=1){
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

  const mxArray * adjacency_mx = prhs[0];
  
  const int * size_A = mxGetDimensions(adjacency_mx);
  int n_node = size_A[0];
  double * adjacency_0 = mxGetPr(adjacency_mx);
  int nnz_adjacency = 0;
  Hash_UInt32_Double adjacency;
  for(int i=0; i<n_node; i++)
    for(int j=0; j<n_node; j++)
      if(i!=j){
        adjacency[get_adj_id_32(i,j)] = fabs(adjacency_0[i + j*n_node]);
        nnz_adjacency++;
      }
      else
        adjacency[get_adj_id_32(i,j)] = 0;
  nnz_adjacency /= 2;
  
  // create a sorted list of adjacencies
  Weighted_Edge * adj_list = (Weighted_Edge *) new
    Weighted_Edge[nnz_adjacency];
  {
    unsigned int i, j, k=0;
    for(i=0; i<n_node; i++){
      for(j=i+1; j<n_node; j++){
        if(adjacency[get_adj_id_32(i,j)]==0)
          continue;
        adj_list[k].w = -adjacency[get_adj_id_32(i,j)];
        adj_list[k].i = i;
        adj_list[k].j = j;
        k++;
      }
    }

    std::sort(adj_list, adj_list+nnz_adjacency);
    for(k=0; k<nnz_adjacency; k++)
      adj_list[k].w *= -1;
  }
  
  // get all 3-cliques in G
  unsigned int n_triplet = 0;
  Hash_UInt64_UInt32 triplet_id;
  n_triplet = 0;
  {
    unsigned int i, j, k;
    unsigned long int l;
    for(i=0; i<n_node; i++){
      for(j=i+1; j<n_node; j++){
        if(adjacency[get_adj_id_32(i,j)]==0)
          continue;
        for(k=0; k<n_node; k++){
          if(k==i || k==j)
            continue;
          if(adjacency[get_adj_id_32(i,k)]==0 ||
             adjacency[get_adj_id_32(j,k)]==0)
            continue;

          l = get_triplet_word(k, i, j);
          if(triplet_id[l]==0){
            n_triplet++;
            triplet_id[l] = n_triplet;
          }
        }
      }
    }
  }
  mexPrintf("n_triplet: %d\n", n_triplet);
  if(n_triplet==0){
    plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    return;
  }

  unsigned int * triplet_0 = (unsigned int*) new unsigned int[n_triplet+1];
  unsigned int * triplet_1 = (unsigned int*) new unsigned int[n_triplet+1];
  unsigned int * triplet_2 = (unsigned int*) new unsigned int[n_triplet+1];
  {
    Hash_UInt64_UInt32::iterator hp;
    unsigned long int i, l, m;
    for(hp=triplet_id.begin(); hp!=triplet_id.end(); hp++){
      i = (*hp).second;
      l = (*hp).first;
      m = l>>32;
      triplet_0[i] = (unsigned int) m;
      m = l>>16;
      m = m & 0x000000000000ffff;
      triplet_1[i] = (unsigned int) m;
      m = l;
      m = m & 0x000000000000ffff;
      triplet_2[i] = (unsigned int) m;
#ifdef __DEBUG__
      mexPrintf("<%d, %d, %d>: %d\n", triplet_0[i], triplet_1[i], triplet_2[i],
                i);
#endif
    }
  }

  // get adjacency matrix for H
  unsigned int n_vertex = triplet_id.size();
  std::vector<Pair_UInt32> edges;
  {
    Hash_UInt64_UInt32::iterator hp;
    unsigned long int i, l, m, n1, n2, n3;
    unsigned int a, b, c, d;
    for(hp=triplet_id.begin(); hp!=triplet_id.end(); hp++){
      i = (*hp).second;
      l = (*hp).first;
      m = l>>32;
      a = (unsigned int) m;
      m = l>>16;
      m = m & 0x000000000000ffff;
      b = (unsigned int) m;
      m = l;
      m = m & 0x000000000000ffff;
      c = (unsigned int) m;
      for(d=0; d<n_node; d++){
        if(d==a || d==b || d==c)
          continue;
        n1 = get_triplet_word(d, a, b);
        n2 = get_triplet_word(d, a, c);
        n3 = get_triplet_word(d, b, c);

        if(adjacency[get_adj_id_32(d,a)]!=0 &&
           adjacency[get_adj_id_32(d,b)]!=0){
          if(n1>l){
            edges.push_back(Pair_UInt32(triplet_id[l], triplet_id[n1]));
          }
        }
        if(adjacency[get_adj_id_32(d,a)]!=0 &&
           adjacency[get_adj_id_32(d,c)]!=0){
          if(n2>l){
            edges.push_back(Pair_UInt32(triplet_id[l], triplet_id[n2]));
          }
        }
        if(adjacency[get_adj_id_32(d,c)]!=0 &&
           adjacency[get_adj_id_32(d,b)]!=0){
          if(n3>l){
            edges.push_back(Pair_UInt32(triplet_id[l], triplet_id[n3]));
          }
        }
      }
    }
  }
  mexPrintf("n_edge: %d\n", edges.size());
  if(edges.size()==0){
    plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    return;
  }

  // compute maximal induced acyclic graph
  Hash_UInt32_UInt32 adj_flag;
  // Simple region growing
  {
    // construct graph
    Graph g(n_triplet);
    {
      graph_traits<Graph>::edge_descriptor e;
      bool is_inserted;
      std::vector<Pair_UInt32>::iterator vp;
      for(vp=edges.begin(); vp!=edges.end(); vp++){
        tie(e, is_inserted) = add_edge((*vp).first, (*vp).second, g);
#ifdef __DEBUG__
        mexPrintf("(%d) -- (%d)\n", (*vp).first, (*vp).second);
#endif
      }
    }

    // greedily add adjacencies
    Hash_UInt32_UInt32 vertex_flag;
    Merge_Sets_H<unsigned int, std::hash<unsigned int>, equi>
      vertex_sets(NULL);
    // first initialize flags
    graph_traits<Graph>::vertex_iterator vs, ve;
    tie(vs, ve) = vertices(g);
    for(; vs!=ve; vs++){
      vertex_flag[*vs] = 0;
      vertex_sets.add_new_set_inc(*vs);
    }
    {
      for(int i=0; i<n_node; i++)
        for(int j=i+1; j<n_node; j++)
          adj_flag[get_adj_id_32(i,j)] = 0;
    }

    // now start adding
    graph_traits<Graph>::out_edge_iterator e_vs, e_vs_end, e_vs1, e_vs1_end;
    
    for(unsigned int adj_list_id = 0; adj_list_id<nnz_adjacency;
        adj_list_id++){
      unsigned int i = adj_list[adj_list_id].i;
      unsigned int j = adj_list[adj_list_id].j;
      if(i==j)
        continue;

      unsigned int a_ij = get_adj_id_32(i,j);
      if(adj_flag[a_ij]>0) // already processed - skip
        continue;
      
#ifdef __DEBUG__
      mexPrintf("i: %d, j: %d\n", i, j);
#endif
      
      // check if a_ij will create a cycle with already included
      // vertices
      //
      // First, get all the neighbors of vertices that will be created
      // by adding this adjacency.
      bool is_present_cycle = false;
      std::set<unsigned int,ltui> neigh_union;
#ifdef __DEBUG__
      mexPrintf("vertices  already created:");
#endif
      for(unsigned int k=0; k<n_node; k++){
        unsigned int a_ik = get_adj_id_32(i, k);
        unsigned int a_jk = get_adj_id_32(j, k);
      
        if(adjacency[a_ik]==0 || adjacency[a_jk]==0)
          continue;
        if(adj_flag[a_ik]==1 && adj_flag[a_jk]==1){
          graph_traits<Graph>::adjacency_iterator v_vs, v_vs_end;
          tie(v_vs, v_vs_end) = adjacent_vertices(
            triplet_id[get_triplet_word(i, j, k)], g);
          for(; v_vs!=v_vs_end; v_vs++)
            if(vertex_flag[*v_vs]==1){
#ifdef __DEBUG__
              mexPrintf("%u, ", *v_vs);
#endif
              if(neigh_union.find(*v_vs)!=neigh_union.end()){
                is_present_cycle = true;
                break;
              }
              neigh_union.insert(*v_vs);
            }
        }
        if(is_present_cycle)
          break;
      }
#ifdef __DEBUG__
      mexPrintf("\n");
#endif
      
      // Next, see if any pairs of neighbors are in the same merge set
      if(!is_present_cycle){
        std::set<unsigned int,ltui>::iterator np;
        Hash_UInt32_UInt32 n_adam;
        for(np=neigh_union.begin(); np!=neigh_union.end(); np++)
          n_adam[vertex_sets.get_adam(*np)]++;
        Hash_UInt32_UInt32::iterator hp;
        for(hp=n_adam.begin(); hp!=n_adam.end(); hp++)
          if((*hp).second>1){
            is_present_cycle=true;
            break;
          }
      }
      
      if(!is_present_cycle){
        adj_flag[a_ij] = 1;
        for(unsigned int k=0; k<n_node; k++){
          unsigned int a_ik = get_adj_id_32(i, k);
          unsigned int a_jk = get_adj_id_32(j, k);
      
          if(adjacency[a_ik]==0 || adjacency[a_jk]==0)
            continue;
          if(adj_flag[a_ik]==1 && adj_flag[a_jk]==1){
            unsigned int v_ijk = triplet_id[get_triplet_word(i, j, k)];
#ifdef __DEBUG__
            mexPrintf("vertex (%d,%d,%d) added\n", i, j, k);
#endif
            vertex_flag[v_ijk] = 1;
            tie(e_vs, e_vs_end) = out_edges(v_ijk, g);
            while(e_vs!=e_vs_end){
              unsigned int n = target(*e_vs, g);
              if(vertex_flag[n]==1)
                vertex_sets.merge(v_ijk,n);
              e_vs++;
            }
          }
        }
      }
      else
        adj_flag[a_ij] = 2;

#ifdef __DEBUG__
      {
        Hash_UInt32_UInt32::iterator hp;
        for(hp=vertex_flag.begin(); hp!=vertex_flag.end(); hp++){
          mexPrintf("vertex_sets[%u]: %u\n", (*hp).first, 
                    vertex_sets.get_adam((*hp).first));
        }
      }
#endif
    }
  }

  
  // output triplet-vertices retained in induced sub-graph
  plhs[0] = mxCreateDoubleMatrix(n_node, n_node, mxREAL);
  {
    Hash_UInt32_Double::iterator hp;
    double * adjacency_new = mxGetPr(plhs[0]);
    for(hp=adjacency.begin(); hp!=adjacency.end(); hp++){
      if((*hp).second!=0 && adj_flag[(*hp).first]==1){
        unsigned int i = (*hp).first >> 16;
        unsigned int j = (*hp).first & 0x0000ffff;
        adjacency_new[i + j*n_node] = 1;
        adjacency_new[j + i*n_node] = 1;
      }
    }

    for(unsigned int i=0; i<n_node; i++)
      adjacency_new[i + i*n_node] = 1;
  }

  delete [] triplet_0;
  delete [] triplet_1;
  delete [] triplet_2;
  
  return;
}
