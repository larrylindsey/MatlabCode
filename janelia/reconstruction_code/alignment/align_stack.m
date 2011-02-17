function align_stack(config)
% Global alignment
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

if(~isfield(config.stack, 'fold') || ...
    ~isfield(config.stack.fold, 'is_considered_in_alignment') || ...
    ~config.stack.fold.is_considered_in_alignment)
  align_stack_SIFT_affine_whole_tile_2(config);
else
  align_stack_dmesh_affine_tile_patches(config)
end

return
end
