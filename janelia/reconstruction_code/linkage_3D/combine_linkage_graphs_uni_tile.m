function combine_linkage_graphs_uni_tile(config)
% combine_linkage_graphs_uni_tile(config)
% For one tile per plane stacks, combine the linkage graphs of consecutive
% section pairs into one linkage graph.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  10222008  init code
%

stack_config = config.stack;
linkage_config = config.linkage;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
apply_config = linkage_config.apply;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Merging the linkage graphs from all sections ..\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
links_3D_t = {};
for layer_id = 1:length(stack_config.case_ids)-1
  case_id = stack_config.case_ids(layer_id);
  fprintf('%d ', case_id);
  if(mod(layer_id,20)==0)
    fprintf('\n');
  end
  load2([get_reconstruction_dir(config), linkage_config.dir, 'links_3D_p', '.', ...
    num2str(case_id), '.', model_config.type, '.', feature_config.type, ...
    apply_config.model_suffix, config.segmentation_choose.choice.seg_suffix, ...
    '.mat'], 'links_3D_p');
  
  links_3D_t{layer_id} = links_3D_p{1,1}; %#ok<USENS>
end

save2([get_region_dir(config), 'links_3D_t', '.', ...
  model_config.type, '.', feature_config.type, apply_config.model_suffix, ...
  config.segmentation_choose.choice.seg_suffix, '.mat'], 'links_3D_t');

fprintf('\ndone.\n');
return
end
