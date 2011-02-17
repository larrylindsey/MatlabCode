function align_stack_tile_pair_in_plane(config)
% align_stack_tile_pair_in_plane(config)
% Compute transformations between pairs of tiles within a section.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

switch(config.align.in_section_tile_pair_align.method)
  case 'SIFT'
    align_stack_SIFT_affine_tile_pair_in_plane(config);
  case 'deformable_mesh'
    if(~config.job.is_stand_alone)
      align_stack_deformable_mesh_tile_pair_in_plane(config);
    else
      align_stack_deformable_mesh_tile_pair_in_plane_SAB(config);
    end
  otherwise
    error('Within plane tile pair alignment method not recognized');
end

return
end
