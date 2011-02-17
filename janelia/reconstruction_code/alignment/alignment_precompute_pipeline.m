function alignment_precompute_pipeline(config, sub_module_id)

precompute_config = config.align.precompute;

% For normalized cross correlation based alignment within sections.
% Transformation constrained to translation.
if(isfield(precompute_config, 'norm_cross_corr_in_plane') && ...
    isfield(precompute_config.norm_cross_corr_in_plane, 'is_enabled') && ...
    precompute_config.norm_cross_corr_in_plane.is_enabled)
  switch(sub_module_id)
    case 1
      align_norm_cross_corr_tile_pair_in_plane(config);
  end  
end

% For normalized cross correlation based alignment in adjacent sections.
% Transformation constrained to translation.
if(isfield(precompute_config, 'norm_cross_corr_inter_plane') && ...
    isfield(precompute_config.norm_cross_corr_inter_plane, 'is_enabled') && ...
    precompute_config.norm_cross_corr_inter_plane.is_enabled)
  switch(sub_module_id)
    case 1
      align_norm_cross_corr_tile_pair_inter_plane(config);
  end  
end

% For SIFT: feature point detection, matching within section and across
% adjacent sections
if(isfield(precompute_config, 'SIFT') && ...
    isfield(precompute_config.SIFT, 'is_enabled') && ...
    precompute_config.SIFT.is_enabled)
  switch(sub_module_id)
    case 1
      align_stack_SIFT_feature_points(config);
    case 2
      align_stack_SIFT_in_plane_matches(config);
    case 3
      align_stack_SIFT_inter_plane_matches(config);
  end
end

% For deformable mesh based computation of correspondences in case SIFT
% fails. This is meant for within plane stitching and ignores folds
if(isfield(precompute_config, 'deformable_mesh_in_plane_overlap') && ...
    isfield(precompute_config.deformable_mesh_in_plane_overlap, 'is_enabled') && ...
    precompute_config.deformable_mesh_in_plane_overlap.is_enabled)
  switch(sub_module_id)
    case 1
      correspond_tile_pair_within_section_deformable_mesh(config);
  end  
end

% For deformable mesh based computation of correspondences in case SIFT
% fails. This is meant for global alignment and considers folds
if(isfield(precompute_config, 'deformable_mesh_in_plane_overlap_with_fold') && ...
    isfield(precompute_config.deformable_mesh_in_plane_overlap_with_fold, 'is_enabled') && ...
    precompute_config.deformable_mesh_in_plane_overlap_with_fold.is_enabled)
  switch(sub_module_id)
    case 1
      if(~config.job.is_stand_alone)
        align_stack_deformable_mesh_tile_pair_in_plane(config);
      else
        align_stack_deformable_mesh_tile_pair_in_plane_SAB(config);
      end
  end  
end

return
end
