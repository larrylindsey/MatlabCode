// Get minimum spanning tree given a directed graph's adjacency matrix
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
//

#include <mex.h>
#include <utility>                       // for std::pair
#include <algorithm>                     // for std::for_each
#include <boost/utility.hpp>             // for boost::tie
#include <boost/graph/graph_traits.hpp>  // for boost::graph_traits
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/properties.hpp>

namespace std{
  using namespace __gnu_cxx;
}

// type of the graph: unsorted sequence of edges, unsorted set of vertices,
// undirected graph, no vertex property and edge property Boundary_Stat
typedef boost::property<boost::edge_weight_t, double> DistanceProperty;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, DistanceProperty> Graph;

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    mexPrintf("input params:\n");
    mexPrintf("\t1. MxM double sparse matrix\n");
    mexPrintf("output:\n");
    mexPrintf("\t1. (M-1)x2 uint32 edge list in MST\n");
    return;
  }
  if(nrhs!=1){
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

  mexPrintf("START: get_minimum_spanning_tree\n");
  const mxArray * A_mx = prhs[0];

  const double * A = mxGetPr(A_mx);
  mwIndex * ir = mxGetIr(A_mx);
  mwIndex * jc = mxGetJc(A_mx);
  const unsigned int n_node = (unsigned int) mxGetN(A_mx);
  mexPrintf("number of nodes: %d\n", n_node);
  
  plhs[0] = mxCreateNumericMatrix(n_node-1, 2, mxUINT32_CLASS,
                                 mxREAL);
  unsigned int * mst_edges = (unsigned int *) mxGetPr(plhs[0]);

  // construct from adjacency matrix
  Graph g(n_node);
  boost::property_map<Graph, boost::edge_weight_t>::type weight = 
    boost::get(boost::edge_weight, g);
  {
    unsigned int n;
    unsigned int n_neigh, edge_id, neigh_id, i, node_id;
    for(n = 0; n<n_node; n++){
      for(edge_id=jc[n]; edge_id<jc[n+1]; edge_id++){
	neigh_id = ir[edge_id];
	bool is_inserted;
	boost::graph_traits<Graph>::edge_descriptor e;
	boost::tie(e, is_inserted) = boost::add_edge(n, neigh_id, g);
	weight[e] = A[edge_id];
#ifdef __DEBUG__
	mexPrintf("w(%d, %d): %g\n", n+1, neigh_id+1, weight[e]);
#endif
      }
    }
  }

  // compute minimum spanning tree
  std::vector<boost::graph_traits<Graph>::edge_descriptor> mste(n_node-1);
  boost::kruskal_minimum_spanning_tree(g, mste.begin());

  // copy result to output
  {
    std::vector<boost::graph_traits<Graph>::edge_descriptor>::iterator e;
    int i=0;
    for(e=mste.begin(); e!=mste.end(); e++, i++){
      mst_edges[i] = source(*e, g)+1; // 0 index in mex
      mst_edges[i+n_node-1] = target(*e, g)+1; // 0 index in mex
    }
  }
  mexPrintf("STOP: get_minimum_spanning_tree\n");
  return;
}
