function display_aligned_tile_patches(config)
% display_aligned_tile_patches(config)
% Display aligned tile images with piecewise affine transformations, e.g.,
% in case there are folds.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  08312008  init. code
% v1  12272008  modified for piecewise affine transformations.
%

stack_config = config.stack;
display_config = config.align.global_align.display;
if(~isfield(display_config, 'is_verbose'))
  display_config.is_verbose = true;
end
if(display_config.is_verbose)
  fprintf('Display globally aligned tile patches\n');
end
switch(config.align.global_align.method)
  case 'SIFT'
    tform_dir = [get_region_dir(config), config.SIFT.dir];
  case 'deformable_mesh'
    tform_dir = [get_region_dir(config), config.deformable_mesh.dir];
end

if(isfield(display_config, 'save_as_stack_file_name') && ...
    ~isempty(display_config.save_as_stack_file_name))
  if(exist(display_config.save_as_stack_file_name, 'file')==2)
    delete(display_config.save_as_stack_file_name);
  end
end

if(display_config.is_verbose)
  fprintf('Putting all patches in global coordinates ...\n');
end
minx = inf;
miny = inf;
maxx = -inf;
maxy = -inf;
tforms = {};
xdata = {};
ydata = {};
images = get_image_from_stack(config, stack_config.case_ids(1));
if(isfield(display_config, 'display_scale') && ~isempty(display_config.display_scale))
  images{1} = imresize(images{1}, 1/display_config.display_scale);
end
size_image = size(images{1});
% for layer_id = [1, 2, length(stack_config.case_ids)-1, length(stack_config.case_ids)]
for layer_id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(layer_id);
  if(display_config.is_verbose)
    fprintf('plane: %d\n', case_id);
  end
  
  image_prefixes  = get_image_prefixes_subdirs(config, case_id);

  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))
      continue;
    end
    if(display_config.is_verbose)
      fprintf('tile: %d\n', tile_id);
    end
    load2([tform_dir, image_prefixes{tile_id}, ...
      '.patchwise_affine_transforms.mat'], 'transform');
    n_patch = size(transform.transform, 2);
    for patch_id = 1:n_patch
      t = reshape(transform.transform(:,patch_id), [2 3]);
      if(isfield(display_config, 'display_scale') && ~isempty(display_config.display_scale))
        t(:,3) = t(:,3)/display_config.display_scale;
      end
      tforms{layer_id, tile_id, patch_id} = ...
        maketform('affine', t'); %#ok<AGROW>
      [junk, xdata{layer_id, tile_id, patch_id}, ydata{layer_id, tile_id, patch_id}] = ...
        imtransform(ones(size_image), tforms{layer_id, tile_id, patch_id}); %#ok<AGROW>
      minx = min([minx, xdata{layer_id, tile_id, patch_id}]);
      miny = min([miny, ydata{layer_id, tile_id, patch_id}]);
      maxx = max([maxx, xdata{layer_id, tile_id, patch_id}]);
      maxy = max([maxy, ydata{layer_id, tile_id, patch_id}]);
    end
  end
end
if(display_config.is_verbose)
  fprintf('done.\n');
end

% for layer_id = [1, 2, length(stack_config.case_ids)-1, length(stack_config.case_ids)]
for layer_id = 1:length(stack_config.case_ids)
  for tile_id = 1:size(xdata,2)
    for patch_id = 1:size(xdata,3)
      if(isempty(xdata{layer_id, tile_id, patch_id}))
        continue;
      end
      xdata{layer_id, tile_id, patch_id} = ...
        round(xdata{layer_id, tile_id, patch_id}-minx+1); %#ok<AGROW>
      ydata{layer_id, tile_id, patch_id} = ...
        round(ydata{layer_id, tile_id, patch_id}-miny+1); %#ok<AGROW>
    end
  end
end
maxx = round(maxx - minx + 1);
maxy = round(maxy - miny + 1);
if(display_config.is_verbose)
  fprintf('done.\n');
end

if(display_config.is_verbose)
  fprintf('Rendering aligned patches ...\n');
end
fold_dir = get_fold_dir(config);
% for layer_id = [1, 2, length(stack_config.case_ids)-1, length(stack_config.case_ids)]
for layer_id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(layer_id);
  if(display_config.is_verbose)
    fprintf('plane: %d\n', case_id);
  end
  [images, image_prefixes]  = get_image_from_stack(config, case_id);

  tile_image = zeros(maxy, maxx, 3);
  t = 0;
  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))
      continue;
    end
    if(display_config.is_verbose)
      fprintf('tile: %d\n', tile_id);
    end
    t = t+1;
    if(isfield(display_config, 'display_scale') && ~isempty(display_config.display_scale))
      images{tile_id} = imresize(images{tile_id}, 1/display_config.display_scale, 'bilinear');
    end
    if(isfield(display_config, 'filter_version') && ~isempty(display_config.filter_version))
      images{tile_id} = filter_image(images{tile_id}, display_config.filter_version);
    end
    if(display_config.is_inverted_display)
      images{tile_id} = 1-images{tile_id};
    end
    load2([tform_dir, image_prefixes{tile_id}, ...
      '.patchwise_affine_transforms.mat'], 'transform');
    file_name_prefix = [fold_dir, image_prefixes{tile_id}, '.fold_mask'];
    fold_mask = load_fold_mask(file_name_prefix, config);
    if(isfield(display_config, 'display_scale') && ~isempty(display_config.display_scale))
      fold_mask = imresize(fold_mask, 1/display_config.display_scale, 'nearest');
    end
    n_patch = size(transform.transform, 2);
    for patch_id = 1:n_patch
      image = images{tile_id} .* (fold_mask==patch_id);
      image_t = imtransform(image, tforms{layer_id, tile_id, patch_id});
      tile_image(...
        ydata{layer_id, tile_id, patch_id}(1):ydata{layer_id, tile_id, patch_id}(2), ...
        xdata{layer_id, tile_id, patch_id}(1):xdata{layer_id, tile_id, patch_id}(2), ...
        mod(t,3)+1) = ...
        max(tile_image(...
        ydata{layer_id, tile_id, patch_id}(1):ydata{layer_id, tile_id, patch_id}(2), ...
        xdata{layer_id, tile_id, patch_id}(1):xdata{layer_id, tile_id, patch_id}(2), ...
        mod(t,3)+1), image_t);
    end
  end
  if(isfield(display_config, 'is_in_grayscale') && ...
      display_config.is_in_grayscale)
    tile_image = max(tile_image, [], 3);
  end
  if(~isfield(display_config, 'save_as_stack_file_name') || ...
      isempty(display_config.save_as_stack_file_name))
    figure; imshow(tile_image);
    title(['plane: ', num2str(case_id)]);
  else
    figure(1); imshow(tile_image);
    title(['plane: ', num2str(case_id)]);
    if(~isfield(display_config, 'display_scale')  || isempty(display_config.display_scale))
      tile_image = imresize(tile_image, 0.25, 'bilinear');
    end
    tile_image(tile_image>1) = 1;
    tile_image(tile_image<0) = 0;
    imwrite(tile_image, sprintf(display_config.save_as_stack_file_name, case_id), ...
      'Compression', 'none', 'WriteMode', 'append');
  end
end

return
end
