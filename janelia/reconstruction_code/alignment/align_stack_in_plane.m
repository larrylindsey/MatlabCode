function align_stack_in_plane(config)
% align_stack_in_plane(config)
% Compute transformations to align tiles within a section
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

in_plane_align_config = config.align.in_section_align;

switch(in_plane_align_config.method)
  case 'SIFT'
    error('Deprecated');
%     align_stack_SIFT_affine_in_plane(config);
  case 'deformable_mesh'
    align_stack_dmesh_affine_in_plane(config);
  otherwise
    error('Within plane alignment method not recognized');
end

return
end
