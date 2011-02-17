function combine_linkage_graphs(config)
% combine_linkage_graphs(config)
% Combine the linkage graphs of consecutive
% section pairs into one linkage graph.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  10222008  init code
%
if(isfield(config.stack, 'image_structure') && ~isempty(config.stack.image_structure))
  combine_linkage_graphs_multi_tile(config);
else
  combine_linkage_graphs_uni_tile(config);
end

return
end
