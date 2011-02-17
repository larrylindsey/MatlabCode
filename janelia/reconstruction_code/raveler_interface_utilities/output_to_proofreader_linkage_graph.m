function seg_2_body_map = output_to_proofreader_linkage_graph(...
  links_3D, superpixel_2_seg_map, label_pairs_overlap, config)
% output_to_proofreader_linkage_graph(links_3D, links_3D_dummy, ...
% superpixel_2_seg_map, config)
% Save linkage graph for proofreader: for MATLAB gui save linkage graph,
% for Raveler save the segment-to-body-map.
% Inputs: 
%   links_3D          linkage graph
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

switch(config.proofreader.method)
  case 'matlab_gui'
    linkage_config = config.linkage;
    model_config = linkage_config.model;
    feature_config = linkage_config.feature;
    apply_config = linkage_config.apply;
    save2([get_to_be_proofread_dir(config), 'links_3D',  '.',  model_config.type, '.', ...
      feature_config.type, apply_config.model_suffix, ...
      config.segmentation_choose.choice.seg_suffix, '.mat'], 'links_3D', ...
      'label_pairs_overlap');
  case 'Raveler'
    seg_2_body_map = output_to_raveler_linkage_graph(...
      links_3D, superpixel_2_seg_map, label_pairs_overlap, config);
  otherwise
    error('Proofreader method not recognized');
end

return
end
