function display_aligned_whole_tile(config)
% display_aligned_whole_tile(config)
% Display aligned tile images
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  08312008  init. code
%

stack_config = config.stack;
if(~isfield(config.align.global_align, 'display'))
  % for backward compatibility
  display_config = config.align.global_align.SIFT;
else
  display_config = config.align.global_align.display;
end
sift_dir = [get_region_dir(config), config.SIFT.dir];

if(isfield(display_config, 'save_as_stack_file_name') && ...
    ~isempty(display_config.save_as_stack_file_name))
  if(exist(display_config.save_as_stack_file_name, 'file')==2)
    delete(display_config.save_as_stack_file_name);
  end
end

minx = inf;
miny = inf;
maxx = -inf;
maxy = -inf;
tforms = {};
xdata = {};
ydata = {};
images = get_image_from_stack(config, stack_config.case_ids(1));
if(~isfield(display_config, 'display_scale') || isempty(display_config.display_scale))
  display_config.display_scale = 1;
end
image_s = imresize(images{1}, 1/display_config.display_scale);
size_image = size(image_s);
for layer_id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(layer_id);
  
  image_prefixes  = get_image_prefixes_subdirs(config, case_id);

  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))
      continue;
    end
    load2([sift_dir, image_prefixes{tile_id}, '.affine_transforms.mat'], 'transform');
    transform(5:6) = transform(5:6)/display_config.display_scale;
    tforms{layer_id, tile_id} = ...
      maketform('affine', reshape(transform, [2 3])'); %#ok<AGROW>
    [junk, xdata{layer_id, tile_id}, ydata{layer_id, tile_id}] = ...
      imtransform(ones(size_image), tforms{layer_id, tile_id}); %#ok<AGROW>
    minx = min([minx, xdata{layer_id, tile_id}]);
    miny = min([miny, ydata{layer_id, tile_id}]);
    maxx = max([maxx, xdata{layer_id, tile_id}]);
    maxy = max([maxy, ydata{layer_id, tile_id}]);
  end
end

for layer_id = 1:size(xdata,1)
  for tile_id = 1:size(xdata,2)
    if(isempty(xdata{layer_id, tile_id}))
      continue;
    end
    xdata{layer_id, tile_id} = round(xdata{layer_id, tile_id}-minx+1); %#ok<AGROW>
    ydata{layer_id, tile_id} = round(ydata{layer_id, tile_id}-miny+1); %#ok<AGROW>
  end
end
maxx = round(maxx - minx + 1);
maxy = round(maxy - miny + 1);

for layer_id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(layer_id);
  fprintf('plane = %d\n', case_id);
  if(isfield(display_config, 'filter_version') && ~isempty(display_config.filter_version))
    get_images_as_is = true;
  else
    get_images_as_is = false;
  end
  [images, image_prefixes]  = get_image_from_stack(config, case_id, get_images_as_is);

  if(display_config.is_inverted_display)
    tile_image = zeros(maxy, maxx, 3);
  else
    tile_image = zeros(maxy, maxx);
  end
  t = 0;
  for tile_id = 1:size(xdata,2)
    if(isempty(xdata{layer_id, tile_id}))
      continue;
    end
    image_s = imresize(images{tile_id}, 1/display_config.display_scale, 'bilinear');
    if(isfield(display_config, 'filter_version') && ~isempty(display_config.filter_version))
      image_s = filter_image(image_s, display_config.filter_version);
    end
    if(display_config.is_inverted_display)
      image_s = 1-image_s;
    end
    load2([sift_dir, image_prefixes{tile_id}, '.affine_transforms.mat'], 'transform');
    transform(5:6) = transform(5:6)/display_config.display_scale;
    tforms{layer_id, tile_id} = ...
      maketform('affine', reshape(transform, [2 3])'); %#ok<AGROW>
    images_t = imtransform(image_s, tforms{layer_id, tile_id});
    t = t+1;
    if(display_config.is_inverted_display)
      tile_image(ydata{layer_id, tile_id}(1):ydata{layer_id, tile_id}(2), ...
        xdata{layer_id, tile_id}(1):xdata{layer_id, tile_id}(2), mod(t,3)+1) = ...
        max(tile_image(ydata{layer_id, tile_id}(1):ydata{layer_id, tile_id}(2), ...
        xdata{layer_id, tile_id}(1):xdata{layer_id, tile_id}(2), mod(t,3)+1), ...
        images_t);
    else
      tile_image(ydata{layer_id, tile_id}(1):ydata{layer_id, tile_id}(2), ...
        xdata{layer_id, tile_id}(1):xdata{layer_id, tile_id}(2)) = ...
        max(tile_image(ydata{layer_id, tile_id}(1):ydata{layer_id, tile_id}(2), ...
        xdata{layer_id, tile_id}(1):xdata{layer_id, tile_id}(2)), ...
        images_t);
    end
  end
  if(~isfield(display_config, 'save_as_stack_file_name') || ...
      isempty(display_config.save_as_stack_file_name))
    figure; imshow(tile_image);
  else
    figure(1); imshow(tile_image);
    imwrite(tile_image, sprintf(display_config.save_as_stack_file_name, case_id), ...
      'Compression', 'lzw', 'WriteMode', 'append');
  end
end

return
end
