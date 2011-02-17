function generate_cat_superpixel_2_seg_map(config)
% generate_cat_superpixel_2_seg_map(config)
% generate cat and superpixel_2_seg_map
%
% cat{z}'s are the superpixel maps being used in the reconstruction
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
% v1  09302008  merging calls for generating cat and superpixel_2_seg_map
%

fprintf('START: generate_cat_superpixel_2_seg_map\n');
if(isfield(config.proofreader.export, 'render_externally') && ...
    config.proofreader.export.render_externally)
  generate_superpixel_2_seg_map_multi_tile_patches(config);
else
  if(isfield(config.stack, 'image_structure') && ~isempty(config.stack.image_structure))
    if(~isfield(config.stack, 'fold') || ...
        ~isfield(config.stack.fold, 'is_considered_in_segmentation') || ...
        ~config.stack.fold.is_considered_in_segmentation)
      generate_cat_superpixel_2_seg_map_multi_tile(config);
    else
      generate_cat_superpixel_2_seg_map_multi_tile_patches(config);
    end
  else
    generate_cat_uni_tile(config);
    generate_superpixel_2_seg_map_uni_tile(config);
  end
end
fprintf('STOP: generate_cat_superpixel_2_seg_map\n');

return
end
