function output_to_proofreader_superpixel_2_seg_map(superpixel_2_seg_map, ...
  superpixel_ids_sets, config)
% output_to_proofreader_superpixel_2_seg_map(superpixel_2_seg_map, ...
%   superpixel_ids_sets, config)
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

switch(config.proofreader.method)
  case 'matlab_gui'
    save2([get_to_be_proofread_dir(config), 'superpixel_2_seg_map', ...
      config.segmentation_choose.choice.seg_suffix, '.mat'], ...
      'superpixel_2_seg_map');
  case 'Raveler'
    output_to_raveler_superpixel_2_seg_map(superpixel_2_seg_map, ...
      superpixel_ids_sets, config);
  otherwise
    error('Proofreader method not recognized');
end

return
end
