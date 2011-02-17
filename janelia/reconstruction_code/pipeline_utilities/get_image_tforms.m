function tforms = get_image_tforms(config, case_id)
% tforms = get_image_tforms(config, case_id)
% Retrieves transforms for images in a section.
%
% Input: 
%   config    config datastructure of the reconstruction
%   case_id   z-plane for which images are to be retrived
% Output:
%   tforms    cell array of affine transforms, one for each tile in a
%               section.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  07162008  init. code
%

tforms = {};
plane_id = find([config.region.region_structure.planes(:).z] == case_id);

% for now just the 0 tilt-plane
for i = 1:length(config.region.region_structure.planes(plane_id).tilt_planes(1).tiles)
  tile = config.region.region_structure.planes(plane_id).tilt_planes(1).tiles(i);
  tforms{i} = tile.transform; %#ok<AGROW>
end
return;
end
