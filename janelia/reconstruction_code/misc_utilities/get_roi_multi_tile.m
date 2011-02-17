function roi = get_roi_multi_tile(config)
% roi = get_roi_multi_tile(config)
% Compute a "maximalish" square overlap region within slices of an aligned
% stack. There are multiple tiles per section (trakEM xml)
% roi.xmin, xmax ymin ymax
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
%

fprintf('Computing ROI for the stack ..\n');

is_debugging = config.DEBUG;

stack_config = config.stack;
sift_config = config.SIFT;
sift_dir = get_sift_dir(config);

image_dir = get_stack_dir(config);

minx = inf;
miny = inf;
maxx = -inf;
maxy = -inf;
tforms = {};
images_t = {};
xdata = {};
ydata = {};
for layer_id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(layer_id);
  
  [images, image_prefixes] = get_image_from_stack(config, case_id);

  for tile_id = 1:length(images)
    if(isempty(images{tile_id}))
      continue;
    end
    load([sift_dir, image_prefixes{tile_id}, '.sift_transforms.mat']);
    tforms{layer_id, tile_id} = maketform('affine', reshape(transform, [2 3])');
    [images_t{layer_id, tile_id}, xdata{layer_id, tile_id}, ydata{layer_id, tile_id}] = ...
      imtransform(ones(size(images{tile_id})), tforms{layer_id, tile_id}, 'nearest');
    minx = min([minx, xdata{layer_id, tile_id}]);
    miny = min([miny, ydata{layer_id, tile_id}]);
    maxx = max([maxx, xdata{layer_id, tile_id}]);
    maxy = max([maxy, ydata{layer_id, tile_id}]);
  end
end

for layer_id = 1:size(images_t,1)
  for tile_id = 1:size(images_t,2)
    if(isempty(images_t{layer_id, tile_id}))
      continue;
    end
    xdata{layer_id, tile_id} = round(xdata{layer_id, tile_id}-minx+1);
    ydata{layer_id, tile_id} = round(ydata{layer_id, tile_id}-miny+1);
  end
end
maxx = round(maxx - minx + 1);
maxy = round(maxy - miny + 1);

slice_mask = ones(maxy, maxx);
for layer_id = 1:size(images_t,1)
  t = 0;
  tile_image = zeros(maxy, maxx);
  for tile_id = 1:size(images_t,2)
    if(isempty(images_t{layer_id, tile_id}))
      continue;
    end
    t = t+1;
    tile_image(ydata{layer_id, tile_id}(1):ydata{layer_id, tile_id}(2), ...
      xdata{layer_id, tile_id}(1):xdata{layer_id, tile_id}(2)) = ...
      max(tile_image(ydata{layer_id, tile_id}(1):ydata{layer_id, tile_id}(2), ...
      xdata{layer_id, tile_id}(1):xdata{layer_id, tile_id}(2)), ...
      images_t{layer_id, tile_id});
  end
  slice_mask = min(slice_mask, tile_image);
end

if(is_debugging)
  figure(1); imshow(slice_mask);
end;

[py, px] = find(slice_mask);

roi.xmin = min(px);
roi.xmax = max(px);
roi.ymin = min(py);
roi.ymax = max(py);
roi.mask = slice_mask(roi.ymin:roi.ymax, roi.xmin:roi.xmax);

if(is_debugging)
  roi_mask = zeros(size(slice_mask));
  roi_mask(roi.ymin:roi.ymax, roi.xmin:roi.xmax) = 1;
  figure(2); imshow(roi_mask);
end;

fprintf('\ndone.\n');

return

end
