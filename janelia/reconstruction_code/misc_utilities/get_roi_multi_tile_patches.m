function roi = get_roi_multi_tile_patches(config)
% roi = get_roi_multi_tile_patches(config)
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

if(isfield(config, 'DEBUG'))
  is_debugging = config.DEBUG;
else
  is_debugging = false;
end

stack_config = config.stack;

load2([get_global_stitching_param_dir(config), 'canvas.mat'], ...
  'canvas_height', 'canvas_width');

images = get_image_from_stack(config, stack_config.case_ids(1));
size_image = size(images{1});
slice_mask = ones(canvas_height, canvas_width);

if(isfield(config.proofreader, 'align_roi') && ...
    isfield(config.proofreader.align_roi, 'is_identity') && ...
    config.proofreader.align_roi.is_identity)
  roi.xmin = 1;
  roi.xmax = canvas_width;
  roi.ymin = 1;
  roi.ymax = canvas_height;
  roi.mask = slice_mask;
  return;
end

error('to do');

for layer_id = 1:size(xdata,1) %#ok<USENS>
  case_id = stack_config.case_ids(layer_id);
  fprintf('plane: %d\n', case_id);

  tile_image = zeros(canvas_height, canvas_width) - 1;
  for tile_id = 1:size(xdata,2)
    for patch_id = 1:size(xdata,3)
      if(isempty(tforms{layer_id, tile_id, patch_id})) %#ok<USENS>
        continue;
      end
      image_t = imtransform(ones(size_image), tforms{layer_id, tile_id, patch_id}, ...
        'nearest', 'FillValue', -1); %#ok<USENS>
      temp_image = zeros(canvas_height, canvas_width)-1;
      temp_image(...
        ydata{layer_id, tile_id, patch_id}(1):ydata{layer_id, tile_id, patch_id}(2), ...
        xdata{layer_id, tile_id, patch_id}(1):xdata{layer_id, tile_id, patch_id}(2)) = ...
        image_t; %#ok<USENS>
      temp_image(tile_image>=0) = tile_image(tile_image>=0);
      tile_image = temp_image;
    end
  end
  % Compute the convex hull of the union tile images as the mask for the
  % section.
  [py, px] = find(tile_image>0);
  hull_i = convhull(px, py);
  hull_x = px(hull_i);
  hull_y = py(hull_i);
  tile_image = poly2mask(hull_x, hull_y, size(tile_image,1), size(tile_image,2));
  
  % Take intersection of the tile images of the various sections.
  slice_mask = min(slice_mask, tile_image);
end
slice_mask(slice_mask<0) = 0;

if(is_debugging)
  figure(1); imshow(slice_mask);
end;

[py, px] = find(slice_mask>0);

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
