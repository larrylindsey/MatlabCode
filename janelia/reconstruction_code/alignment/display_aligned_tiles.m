function display_aligned_tiles(config)
% Global alignment using SIFT affine
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

if(~isfield(config.stack, 'fold') || ...
    ~isfield(config.stack.fold, 'is_considered_in_alignment') || ...
    ~config.stack.fold.is_considered_in_alignment)
  display_aligned_whole_tile(config)
else
  display_aligned_tile_patches(config)
end

return
end
