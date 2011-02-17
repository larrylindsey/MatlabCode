function output_to_proofreader_stitched_superpixel_maps(cat, config)
% output_to_proofreader_stitched_superpixel_maps(cat, config)
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

fprintf('START: output_to_proofreader_stitched_superpixel_maps\n');

proofread_config = config.proofreader;
export_config = proofread_config.export;
if(isfield(export_config, 'superpixel_map_filter_version') && ...
    ~isempty(export_config.superpixel_map_filter_version))
  fprintf('Applying final filters on superpixel map ...\n');
  for i = 1:legth(cat)
    fprintf('layer id: %d\n', i);
    cat{i} = filter_image(cat{i}, export_config.superpixel_map_filter_version);
  end
  fprintf('done.\n');
end

switch(proofread_config.method)
  case 'matlab_gui'
    save2([get_to_be_proofread_dir(config), 'cat', ...
      config.superpixel_choose.choice.seg_suffix, '.mat'], 'cat', '-v7.3');
  case 'Raveler'
    output_to_raveler_stitched_superpixel_maps(cat, config)
  otherwise
    error('Proofreader method not recognized');
end
fprintf('STOP: output_to_proofreader_stitched_superpixel_maps\n');

return
end
