function sizes = get_image_sizes(config, case_id)
% sizes = get_image_sizes(config, case_id)
% Retrieves sizes of images in a section.
%
% Input: 
%   config    config datastructure of the reconstruction
%   case_id   z-plane for which images are to be retrived
% Output:
%   sizes     cell array of [width, height], one for each tile in a
%               section.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  07162008  init. code
%

sizes = {};
plane_id = find([config.region.region_structure.planes(:).z] == case_id);

% for now just the 0 tilt-plane
for i = 1:length(config.region.region_structure.planes(plane_id).tilt_planes(1).tiles)
  tile = config.region.region_structure.planes(plane_id).tilt_planes(1).tiles(i);
  sizes{i} = [tile.width, tile.height]; %#ok<AGROW>
end
return;
end
