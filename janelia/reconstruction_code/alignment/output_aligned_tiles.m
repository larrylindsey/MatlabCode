function output_aligned_tiles(config)
% output_aligned_tiles(config)
% Store aligned tile images in a global variable "aligned_images"
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  08312008  init. code
%

global aligned_images;

stack_config = config.stack;
sift_dir = [get_region_dir(config), config.SIFT.dir];

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
    load2([sift_dir, image_prefixes{tile_id}, '.sift_transforms.mat'], 'transform');
    tforms{layer_id, tile_id} = ...
      maketform('affine', reshape(transform, [2 3])');
    [images_t{layer_id, tile_id}, xdata{layer_id, tile_id}, ydata{layer_id, tile_id}] = ...
      imtransform(images{tile_id}, tforms{layer_id, tile_id});
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

aligned_images = {};
for layer_id = 1:size(images_t,1)
  for tile_id = 1:size(images_t,2)
    if(isempty(images_t{layer_id, tile_id}))
      continue;
    end
    aligned_images{layer_id, tile_id} = zeros(maxy, maxx);
    aligned_images{layer_id, tile_id}(...
      ydata{layer_id, tile_id}(1):ydata{layer_id, tile_id}(2), ...
      xdata{layer_id, tile_id}(1):xdata{layer_id, tile_id}(2)) = ...
      images_t{layer_id, tile_id};
    figure; imshow(aligned_images{layer_id, tile_id});
  end
end

return
end
