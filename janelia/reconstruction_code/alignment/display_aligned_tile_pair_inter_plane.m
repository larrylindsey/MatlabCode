function display_aligned_tile_pair_inter_plane(config)
% display_aligned_tiles_in_plane(config)
% Display aligned tile images using transformations computed for each
% plane/section separately. This is part of two stage alignment
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  08312008  init. code
% v1  12152008  modified for multiple alignment methods
%

switch(config.align.linkage_align.method)
  case 'SIFT'
    display_aligned_SIFT_affine_tile_pair_inter_plane(config);
  case 'deformable_mesh'
    display_aligned_tile_pair_inter_plane_dmesh(config);
  otherwise
    error('Display inter plane tile pair alignment method not recognized');
end

return
end
