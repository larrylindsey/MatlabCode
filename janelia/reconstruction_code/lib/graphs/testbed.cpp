// testbed for graph utilities

#include <iostream>

#include <digraph.h>

using namespace std;

void print_nodes(vector<gut::Node> nv){
  if(nv.empty())
    return;
  cout << *(nv.begin());
  for(vector<gut::Node>::iterator ni=nv.begin()+1;
      ni!=nv.end(); ni++)
    cout << ',' << (*ni);
}

int main(int argc, char * argv[]){
  gut::Graph G;
  G.add_edge(1,2,1);
  cout << "Is edge 1,2: " << G.is_connected(1,2) << '\n';
  cout << "Is edge 1,2: " << G.is_connected(1,2) << '\n';
  cout << "Is edge 2,1: " << G.is_connected(2,1) << '\n';
  cout << "Is edge 2,1: " << G.is_connected(2,1) << '\n';
  cout << "Is edge 1,2: " << G.get_weight(1,2) << '\n';
  cout << "Is edge 2,1: " << G.get_weight(2,1) << '\n';
  G.add_edge(2,10,13);
  cout << "Is edge 1,2: " << G.is_connected(1,2) << '\n';
  cout << "Is edge 1,2: " << G.is_connected(1,2) << '\n';
  cout << "Is edge 2,1: " << G.is_connected(2,1) << '\n';
  cout << "Is edge 2,1: " << G.is_connected(2,1) << '\n';
  cout << "Is edge 1,2: " << G.get_weight(1,2) << '\n';
  cout << "Is edge 2,1: " << G.get_weight(2,1) << '\n';

  cout << "Out-neighbors of 1: ";
  print_nodes(G.get_neighbors_out(1));
  cout << '\n';
  cout << "In-neighbors of 1: ";
  print_nodes(G.get_neighbors_in(1));
  cout << '\n';
  cout << "In-neighbors of 2: ";
  print_nodes(G.get_neighbors_in(2));
  cout << '\n';
  cout << "In-neighbors of 3: ";
  print_nodes(G.get_neighbors_in(3));
  cout << '\n';
  cout << "In-neighbors of 10: ";
  print_nodes(G.get_neighbors_in(10));
  cout << '\n';

  cout << "Nodes: ";
  print_nodes(G.get_nodes());
  cout << '\n';
  
  return 1;
}
