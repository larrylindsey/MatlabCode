function align_stack_tile_pair_inter_plane(config)
% align_stack_tile_pair_inter_plane(config)
% Compute transformations between pairs of tiles of consequtive
% sections.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

switch(config.align.linkage_align.method)
  case 'SIFT'
    align_stack_SIFT_affine_tile_pair_inter_plane(config);
  case 'deformable_mesh'
    if(~config.job.is_stand_alone)
      align_stack_deformable_mesh_tile_pair_inter_plane(config);
    else
      align_stack_deformable_mesh_tile_pair_inter_plane_SAB(config);
    end
  otherwise
    error('Inter plane tile pair alignment method not recognized');
end

return
end
