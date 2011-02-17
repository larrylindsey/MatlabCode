function compute_filtered_images(config)
% compute_filtered_images(config)
% Compute filtered images
%
%
% Input:
%   config    config datastructure of the reconstruction
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%

filtered_image_config = config.precompute.filtered_image;
if(~isfield(filtered_image_config, 'is_verbose'))
  filtered_image_config.is_verbose = true;
end
if(filtered_image_config.is_verbose)
  fprintf('START: compute_filtered_images\n');
end

stack_config = config.stack;
filtered_image_dir = get_filtered_image_dir(config);

for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  if(filtered_image_config.is_verbose)
    fprintf('case_id: %d\n', case_id);
  end
  [images, image_prefixes, image_sub_dirs, is_to_be_processed] = ...
    get_image_from_stack(config, case_id, true);
  
  for j = 1:length(images)
    image = images{j};
    image_prefix = image_prefixes{j};
    image_sub_dir = image_sub_dirs{j};
    if(filtered_image_config.is_verbose)
      fprintf('%s\n', image_prefix);
    end
    if(is_to_be_processed(j)==0)
      if(filtered_image_config.is_verbose)
        fprintf('is_to_be_processed=false, skipping\n');
      end
      continue;
    end
    
    if(isfield(filtered_image_config, 'filter_version') && ~isempty(filtered_image_config.filter_version))
      image = filter_image2(image, filtered_image_config.filter_version, ...
        config, image_prefix);
    end
    
    check_for_dir([filtered_image_dir, image_sub_dir]);
    imwrite(image, [filtered_image_dir, image_prefix, '.', ...
      regexprep(filtered_image_config.filter_version, '(?<char>[^\\])\', '$<char>'), '.tif']);
  end
end

if(filtered_image_config.is_verbose)
  fprintf('STOP: compute_filtered_images\n');
end

return
end
