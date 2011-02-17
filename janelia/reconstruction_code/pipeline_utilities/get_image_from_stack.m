function [images, image_prefixes, image_sub_dirs, is_to_be_processed, size_images] = ...
  get_image_from_stack(config, case_id, get_images_as_is)
% [images, image_prefixes, image_sub_dirs]= get_image_from_stack(config, case_id)
% Retrieves an image from the stack either from the reconstruction's
% "config" datastructure or from Leginon-format .xml file.
%
% Input: 
%   config    config datastructure of the reconstruction
%   case_id   z-plane for which images are to be retrived
%   get_images_as_is    give raw images without any processing (false)
% Output:
%   images              cell array of retrieved images
%   image_prefixes      cell of image prefix names for saving files
%   image_sub_dirs      cell of subdirectories within which the images are
%                         located. This enables creation of a equivalent
%                         subdirectory in the reconstruction directory.
%   is_to_be_processed  whether an image is to be processed. Useful for 
%                         reconstructing of stack part of which has already
%                         been reconstructed. Modules working on sets of
%                         images will process the set if at least one of
%                         the images has is_to_be_processed set.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  07162008  init. code
%

if(nargin<3)
  get_images_as_is = false;
end

stack_config = config.stack;
image_dir = get_stack_dir(config);

images = {};
image_prefixes = {};
image_sub_dirs = {};
is_to_be_processed = [];
size_images = {};

plane_id = find([config.region.region_structure.planes(:).z] == case_id);

% for now just the 0 tilt-plane
for i = 1:length(config.region.region_structure.planes(plane_id).tilt_planes(1).tiles)
  tile = config.region.region_structure.planes(plane_id).tilt_planes(1).tiles(i);
  images{i} = imread([image_dir, tile.name, stack_config.image_suffix]); %#ok<AGROW>
  if(~get_images_as_is)
    images{i} = im2double(images{i}); %#ok<AGROW>
    if(size(images{i},3)>1)
      images{i} = rgb2gray(images{i}); %#ok<AGROW>
    end
    if(max(images{i}(:))>1)
      images{i} = images{i}/max(images{i}(:)); %#ok<AGROW>
    end
  end
  image_prefixes{i} = tile.name; %#ok<AGROW>
  image_sub_dirs{i} = tile.sub_dir; %#ok<AGROW>
  is_to_be_processed(i) = tile.is_to_be_processed; %#ok<AGROW>
  size_images{i} = [tile.height, tile.width]; %#ok<AGROW>
end
return;
end
