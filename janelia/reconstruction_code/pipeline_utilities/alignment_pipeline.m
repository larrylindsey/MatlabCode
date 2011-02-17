function alignment_pipeline(config, main_module_id)

if(config.is_verbose)
  fprintf('START: alignment_pipeline\n');
end
switch(main_module_id)
  case 120
    align_stack_SIFT_feature_points(config);

  case 130
    align_norm_cross_corr_tile_pair_in_plane(config);
  case 131
    align_norm_cross_corr_tile_pair_inter_plane(config);

  case 150
    align_stack_SIFT_in_plane_matches(config);
  case 151
    correspond_tile_pair_within_section_deformable_mesh(config);
  case 155
    align_stack_tile_pair_in_plane(config);
  case 158
    display_aligned_tile_pair_in_plane(config); % (TODO)

  case 160
    align_stack_SIFT_inter_plane_matches(config);
  case 165
    align_stack_tile_pair_inter_plane(config);
  case 168
    display_aligned_tile_pair_inter_plane(config);

  case 175
    align_stack_in_plane(config);
  case 178
    display_aligned_tiles_in_plane(config);

  case 185
    align_stack(config);
  case 188
    display_aligned_tiles(config);
  otherwise
    error('Alignment: module id not recongized');
end

if(config.is_verbose)
  fprintf('STOP: alignment_pipeline\n');
end
return
end
