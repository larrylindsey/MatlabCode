function roi = get_roi(config)
% roi = get_roi(config)
% Compute a "maximalish" square overlap region within slices of an aligned stack.
% roi.xmin, xmax ymin ymax
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
    roi = get_roi_multi_tile(config);
  else
    roi = get_roi_multi_tile_patches(config);
  end
else
  roi = get_roi_uni_tile(config);
end

return

end
