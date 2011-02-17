function links_3D_p = include_proofread_linkage_graph(...
  image_prefix_1, image_prefix_2, linkage_config, links_3D_p, config)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

file_name_prefix = get_file_name_from_tuple(...
  [get_reconstruction_dir(config), linkage_config.dir], ...
  image_prefix_1, image_prefix_2, 'lkp.');
file_name = [file_name_prefix, '.proofread', '.', ...
  linkage_config.include_proofread.proofread_name, ...
  linkage_config.include_proofread.model_suffix, ...
  config.segmentation_choose.choice.seg_suffix, '.mat'];
linkage_graph_proofread = [];
try
  fprintf('Attempting to load %s ...\n', file_name);
  linkage_graph_proofread = load2(file_name, 'links_3D_p');
catch %#ok<CTCH>
  fprintf('failed.\nProceeding nevertheless.\n');
end
fprintf('successfully loaded.\n');
if(~isempty(linkage_graph_proofread))
  % Remove all references to proofread segments in links_3D_p and
  % add linkage_graph_proofread to it.
  if(~isempty(links_3D_p))
    to_retain = ~ismember(links_3D_p(:,1), ...
      linkage_graph_proofread.links_3D_p(:,1));
    links_3D_p = links_3D_p(to_retain, :);
  end
  if(~isempty(links_3D_p))
    to_retain = ~ismember(links_3D_p(:,2), ...
      linkage_graph_proofread.links_3D_p(:,2));
    links_3D_p = links_3D_p(to_retain, :);
  end
  fprintf('Added proofread linkage graph to links_3D_p.\n');
  links_3D_p = [links_3D_p; linkage_graph_proofread.links_3D_p];
else
  fprintf('However, linkage graph was empty - skipping.\n');
end

return
end
