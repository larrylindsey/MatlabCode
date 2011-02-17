function display_aligned_tiles_affine_in_plane(config)
% display_aligned_tiles_SIFT_affine_in_plane(config)
% Display aligned tile images using transformations computed for each
% plane/section separately. This is part of two stage alignment
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  08312008  init. code
%

stack_config = config.stack;
align_config = config.align.in_section_align;
if(~isfield(align_config, 'is_verbose'))
  align_config.is_verbose = true;
end

if(strcmp(align_config.method, 'SIFT')==1)
  tform_dir = get_sift_dir(config);
  tform_suffix = '.sift_transform.in_plane';
  tform_config = config.align.in_section_align.SIFT;
else
  tform_dir = get_deformable_mesh_dir(config);
  tform_suffix = '.dmesh_affine.in_plane';
  tform_config = config.align.in_section_align.deformable_mesh;
end

if(isfield(config.align.in_section_align, 'display') && ...
    ~isempty(config.align.in_section_align.display))
  display_config = config.align.in_section_align.display;
else
  display_config = tform_config; % for backward compatibility
end

for layer_id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(layer_id);
  if(align_config.is_verbose)
    fprintf('plane:%d\n', case_id);
  end
  
  if(~isfield(display_config, 'filter_version') || ...
      isempty(display_config.filter_version))
    get_images_as_is = false;
  else
    get_images_as_is = true;
  end
  [images, image_prefixes] = get_image_from_stack(config, case_id, get_images_as_is);

  minx = inf;
  miny = inf;
  maxx = -inf;
  maxy = -inf;
  tforms = {};
  images_t = {};
  xdata = {};
  ydata = {};
  if(align_config.is_verbose)
    fprintf('Computing transforms ...\n');
  end
  % Generate a substring for the file name that accounts for the images
  % involved in the joint alignment
  image_set_string = '.';
  for tile_id = 1:length(image_prefixes)
    image_set_string = [image_set_string, image_prefixes{tile_id}]; %#ok<AGROW>
  end
  image_set_string = strrep(image_set_string, '/', '_');
  image_set_string = strrep(image_set_string, '\', '_');
  for tile_id = 1:length(images)
    if(isempty(images{tile_id}))
      continue;
    end
    if(align_config.is_verbose)
      fprintf('tile:%d\n', tile_id);
    end
    load2([tform_dir, image_prefixes{tile_id}, image_set_string, ...
      tform_suffix, '.mat'], 'transform_p');
    if(isfield(display_config, 'display_scale') && ...
        ~isempty(display_config.display_scale))
      images{tile_id} = imresize(images{tile_id}, 1/display_config.display_scale);
      transform_p(5:6) = transform_p(5:6) / display_config.display_scale;
    end
    if(isfield(display_config, 'is_inverted_display') && ...
        display_config.is_inverted_display)
      images{tile_id} = 1-images{tile_id};
    end
    if(isfield(display_config, 'filter_version') && ...
        ~isempty(display_config.filter_version))
      images{tile_id} = filter_image(images{tile_id}, display_config.filter_version);
    end
    tforms{layer_id, tile_id} = ...
      maketform('affine', reshape(transform_p, [2 3])'); %#ok<AGROW>
    [images_t{layer_id, tile_id}, xdata{layer_id, tile_id}, ydata{layer_id, tile_id}] = ...
      imtransform(images{tile_id}, tforms{layer_id, tile_id}); %#ok<AGROW>
    minx = min([minx, xdata{layer_id, tile_id}]);
    miny = min([miny, ydata{layer_id, tile_id}]);
    maxx = max([maxx, xdata{layer_id, tile_id}]);
    maxy = max([maxy, ydata{layer_id, tile_id}]);
  end
  if(align_config.is_verbose)
    fprintf('done\n');
  end

  for tile_id = 1:size(images_t,2)
    if(isempty(images_t{layer_id, tile_id}))
      continue;
    end
    xdata{layer_id, tile_id} = round(xdata{layer_id, tile_id}-minx+1); %#ok<AGROW>
    ydata{layer_id, tile_id} = round(ydata{layer_id, tile_id}-miny+1); %#ok<AGROW>
  end
  maxx = round(maxx - minx + 1);
  maxy = round(maxy - miny + 1);

  if(align_config.is_verbose)
    fprintf('Rendering ... \n');
  end
  if(~isfield(display_config, 'multi_plane') || display_config.multi_plane)
    tile_image = zeros(maxy, maxx, 3);
  else
    tile_image = zeros(maxy, maxx);
  end
  t = 0;
  for tile_id = 1:size(images_t,2)
    if(isempty(images_t{layer_id, tile_id}))
      continue;
    end
    if(~isfield(display_config, 'multi_plane') || display_config.multi_plane)
      t = t+1;
    end
    tile_image(ydata{layer_id, tile_id}(1):ydata{layer_id, tile_id}(2), ...
      xdata{layer_id, tile_id}(1):xdata{layer_id, tile_id}(2), mod(t,3)+1) = ...
      max(tile_image(ydata{layer_id, tile_id}(1):ydata{layer_id, tile_id}(2), ...
      xdata{layer_id, tile_id}(1):xdata{layer_id, tile_id}(2), mod(t,3)+1), ...
      images_t{layer_id, tile_id});
  end
  if(~isfield(display_config, 'save_as_tif_prefix') || ...
      isempty(display_config.save_as_tif_prefix))
    figure; imshow(tile_image);
    title(['plane: ', num2str(case_id)]);
  else
    if(config.is_verbose_figures)
      figure(1); imshow(tile_image);
      title(['plane: ', num2str(case_id)]);
    else
      warning('config.is_verbose_figures is set to false, so not displaying figures.\n'); %#ok<WNTAG>
    end
    if(~isfield(display_config, 'display_scale')  || isempty(display_config.display_scale))
      tile_image = imresize(tile_image, 0.25, 'bilinear');
    end
    imwrite(tile_image, sprintf(display_config.save_as_tif_prefix, case_id), ...
      'Compression', 'none');
  end
  if(align_config.is_verbose)
    fprintf('done\n');
  end
end

return
end
