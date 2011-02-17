function generate_coarse_al_without_roi(config)
% generate_coarse_al_without_roi(config)
% generate coarse al without proofreader roi for browsing through the stack
% inorder to define the proofreader roi. Use before generating the final al
% and superpixel mat.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
%

if(isfield(config.stack, 'image_structure') && ~isempty(config.stack.image_structure))
  if(~isfield(config.stack, 'fold') || ...
      ~isfield(config.stack.fold, 'is_considered_in_segmentation') || ...
      ~config.stack.fold.is_considered_in_segmentation)
    generate_al_multi_tile_coarse(config);
  else
    generate_al_multi_tile_patches_coarse(config);
  end
else
  generate_al_uni_tile_coarse(config);
end

return
end
