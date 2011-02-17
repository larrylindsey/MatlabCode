function component_ids = get_connected_components(A, threshold)
% get_connected_components(A, threshold)
% Compute connected components in an undirected weighted graph given its
% weighted adjacency matrix, A, and a threshold. All edges with weight less
% than the threshold are ignored.
%
% Algorithm: simple breadth-first parse of the graph.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

n_node = size(A,1);

% force A to be symmetric
A = sparse(max(A, A') > threshold);

component_ids = zeros(n_node, 1);
processed_flags = zeros(1, n_node);

for adam = 1:n_node
  A(adam,adam)=0;
  if(processed_flags(adam)==1)
    continue;
  end
  processed_flags(adam) = 1;
  component_ids(adam) = adam;
  
  neighbor_ids = find(A(adam,:)==1 & processed_flags==0);
  component_ids(neighbor_ids) = adam;
  node_list = neighbor_ids;
  tail = length(node_list);
  while(tail>0)
    node_id = node_list(tail);
    node_list = node_list(1:end-1);
    tail = tail - 1;
    if(processed_flags(node_id)==1)
      continue;
    end
    A(node_id, node_id)=0;
    
    neighbor_ids = find(A(node_id,:)==1 & processed_flags==0);
    component_ids(neighbor_ids) = adam;
    node_list = [node_list, neighbor_ids];
    tail = length(node_list);
    processed_flags(node_id) = 1;
  end 
end

return
end
