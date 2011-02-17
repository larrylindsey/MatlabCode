// Graph utilities
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI.
//
// v0  03202009   init. code
//

// The adjacency matrix is maintained through hashing, specifically
// stl map utility.

#ifndef __GRAPH_UTILITIES__
#define __GRAPH_UTILITIES__

#include<vector>
#include<set>
#include<map>

namespace gut
{
  typedef unsigned int Node;
  typedef long unsigned int Edge;
  typedef int Weight;

  struct ltlu{
    bool operator()(const long unsigned int l1, const long unsigned int l2) const
      {
        return l1 < l2;
      }
  };
  struct ltu{
    bool operator()(const unsigned int n1, const unsigned int n2) const
      {
        return n1 < n2;
      }
  };

  typedef std::set<Node> Node_Set;
  typedef std::map<const Edge, Weight, ltlu> Edge_Set;

  class Graph{
  public:
    Graph(void);
    ~Graph();

    // Add a node
    void add_node(Node n);

    // Get vector list of nodes
    std::vector<Node> get_nodes(void);
    
    // Add a directed edge from node n1 to node n2
    void add_edge(Node n1, Node n2, Weight w);

    // Find out if an edge exists from node n1 to node n2
    bool is_connected(Node n1, Node n2);

    // Get the weight of edge from node n1 to node n2.
    // Returns 0 if no edge exists
    Weight get_weight(Node n1, Node n2);

    // Get nodes having edges from n
    std::vector<Node> get_neighbors_out(Node n);

    // Get nodes having edges to n
    std::vector<Node> get_neighbors_in(Node n);
    
  private:
    // Set of nodes
    Node_Set N;
    // adjacency matrix and it's reverse for efficiency
    Edge_Set A, A_r;
    Edge get_edge_key(Node n1, Node n2);
  };
}
#endif // __GRAPH_UTILITIES__
