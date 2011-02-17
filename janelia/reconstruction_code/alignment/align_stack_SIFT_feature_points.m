function [frames, descriptors] = align_stack_SIFT_feature_points(config)
% [frames, descriptors] = align_stack_SIFT_feature_points(config)
% Constructs SIFT landmarks for image from the stack either
% from the reconstruction's "config" datastructure or from Leginon-format
% .xml file.
%
% THIS USES getSIFT* LIBRARY
%
% Input:
%   config    config datastructure of the reconstruction
% Output:
%   frames              cell array of SIFT frames (feature point
%                         coordinates)
%   descriptors         cell array of descriptors
%
% Yuriy Mishchenko
% Janelia Farm Research Campus, HHMI.
%
% v0  08252008  init. code
% v1  08312008  modifications for pipeline Shiv N. Vitaladevuni
% v2  09192008  split into SIFT feature point, matching and transform.
%

stack_config = config.stack;
sift_config = config.align.precompute.SIFT.feature_point;

if(~isfield(sift_config, 'is_verbose'))
  sift_config.is_verbose = true;
end
if(~isfield(sift_config, 'is_verbose_figures'))
  sift_config.is_verbose_figures = false;
end
if(sift_config.is_verbose)
  fprintf('START: align_stack_SIFT_feature_points\n');
end

sift_dir = get_sift_dir(config);

for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  fprintf('case_id: %d\n', case_id);
  
  [images, image_prefixes, image_sub_dirs] = get_image_from_stack(...
    config, case_id, true);
  
  if(isfield(sift_config, 'filter_version') && ~isempty(sift_config.filter_version))
    for j = 1:length(images)
      images{j} = filter_image2(images{j}, sift_config.filter_version);
    end
  end
  
  %Generate SIFT frames
  [frames,descriptors]=getSIFTframes(images, sift_config);
  if(sift_config.is_verbose)
    fprintf('done.\n');
  end
  for j=1:size(frames,2)
    if(isempty(image_prefixes{j}))
      continue;
    end
    file_name = [sift_dir, image_prefixes{1,j}, '.sift_landmarks.mat'];
    delete(get_storage_file_name(file_name));
    
    if(isempty(frames{1,j}))
      warning('No SIFT features found for %s', image_prefixes{1,j}); %#ok<WNTAG>
      continue;
    end
    if(sift_config.is_verbose_figures)
      figure;
      imshow(images{1,j});
      hold on;
      plot(frames{1,j}(1,:), frames{1,j}(2,:), '*r');
      hold off;
    end
    check_for_dir([sift_dir, image_sub_dirs{1,j}]);
    frame=frames{1,j}; descriptor=descriptors{1,j}; %#ok<NASGU>
    save2(file_name,'frame','descriptor');
  end
end

if(sift_config.is_verbose)
  fprintf('STOP: align_stack_SIFT_feature_points\n');
end
return
end
