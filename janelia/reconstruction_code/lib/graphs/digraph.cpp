// Graph utilities
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI.
//
// v0  03202009   init. code
//

// The adjacency matrix is maintained through hashing, specifically
// stl map utility.

#include <iostream>
#include <vector>
#include <map>

#include <digraph.h>

namespace gut
{
  Graph::Graph(void){
    return;
  }

  Graph::~Graph(){
    return;
  }
  
  Edge Graph::get_edge_key(Node n1, Node n2){
    Edge l;
    l = n1;
    l <<= 32;
    Edge e;
    e = n2;
    e = e | l;
    return e;
  }

  void Graph::add_node(Node n){
    N.insert(n);
  }

  std::vector<Node> Graph::get_nodes(void){
    std::vector<Node> nv;
    for(Node_Set::iterator ni=N.begin(); ni!=N.end(); ni++)
      nv.push_back(*ni);
    return nv;
  }
  
  void Graph::add_edge(Node n1, Node n2, Weight w){
    N.insert(n1);
    N.insert(n2);
    A[get_edge_key(n1,n2)] = w;
    A_r[get_edge_key(n2,n1)] = w;
  }

  bool Graph::is_connected(Node n1, Node n2){
    return A.find(get_edge_key(n1,n2))!=A.end();
  }

  Weight Graph::get_weight(Node n1, Node n2){
    Edge k = get_edge_key(n1,n2);
    if(A.find(k)==A.end())
      return 0;
    else
      return A[k];
  }

  std::vector<Node> Graph::get_neighbors_out(Node n){
    std::vector<Node> neigh;
    Edge_Set::iterator start = A.lower_bound(get_edge_key(n,0));
    if(start==A.end()){
      return neigh;
    }
    Edge_Set::iterator end = A.upper_bound(get_edge_key(n+1,0));
    Edge_Set::iterator e;

    int i=0;
    Edge l;
    for(e=start; e!=end; e++, i++){
      l = (*e).first;
      l &= 0x00000000ffffffff;
      neigh.push_back((Node) l);
    }
    return neigh;
  }

  std::vector<Node> Graph::get_neighbors_in(Node n){
    std::vector<Node> neigh;
    Edge_Set::iterator start = A_r.lower_bound(get_edge_key(n,0));
    if(start==A_r.end()){
      return neigh;
    }
    Edge_Set::iterator end = A_r.upper_bound(get_edge_key(n+1,0));
    Edge_Set::iterator e;

    int i=0;
    Edge l;
    for(e=start; e!=end; e++, i++){
      l = (*e).first;
      l &= 0x00000000ffffffff;
      neigh.push_back((Node) l);
    }
    return neigh;
  }
}

