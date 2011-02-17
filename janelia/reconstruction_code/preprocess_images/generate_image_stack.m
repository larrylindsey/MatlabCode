function generate_image_stack(config)
% generate_image_stack(config)
% generate image sub-stack TIFFs after filtering.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  03112009  init code
%

stack_config = config.stack;
image_stack_config = config.image_stack;
if(~isfield(image_stack_config, 'is_verbose'))
  image_stack_config.is_verbose = true;
end

if(image_stack_config.is_verbose)
  fprintf('START: generate_image_stack\n');
end
save_file_name = [get_image_stack_dir(config), 'image_stack', ...
  get_image_substack_prefix(config), '.tif'];
if(image_stack_config.is_verbose)
  fprintf('Sub-stack name:%s\n', save_file_name);
end
if(exist(save_file_name, 'file')==2)
  warning('GENERATE_SUBSTACK:OVERWRITE', ...
    ['Overwriting a pre-existing file ', save_file_name]);
  delete(save_file_name);
end
for i = 1:length(stack_config.case_ids)
  plane = stack_config.case_ids(i);
  [images, image_prefixes] = get_image_from_stack(config, plane);
  image = images{1};
  if(image_stack_config.is_verbose)
    fprintf('%s\n', image_prefixes{1});
  end
  
  if(isfield(stack_config, 'roi') && ~isempty(stack_config.roi))
    image = image(stack_config.roi.ymin:stack_config.roi.ymax, ...
      stack_config.roi.xmin:stack_config.roi.xmax, :);
  end
  if(isfield(image_stack_config, 'filter_version') && ...
      ~isempty(image_stack_config.filter_version))
    image = filter_image(image, image_stack_config.filter_version);
  end
  imwrite(image, save_file_name, 'WriteMode', 'append', 'Compression', 'none');
end
if(image_stack_config.is_verbose)
  fprintf('STOP: generate_image_stack.\n');
end
return
end
