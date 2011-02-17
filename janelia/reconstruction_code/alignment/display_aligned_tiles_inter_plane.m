function display_aligned_tiles_inter_plane(config)
% display_aligned_tiles_inter_plane(config)
% Display aligned tile images after two stage alignment. First tiles within
% each section and then pairs of adjacent sections
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  08312008  init. code
%

stack_config = config.stack;
sift_config = config.SIFT;
sift_dir = get_sift_dir(config);

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

  file_name=[sift_dir, 'section.', num2str(case_id), '.sift_transform.global.mat'];
  load(file_name,'transform_g');
  
  for tile_id = 1:length(images)
    if(isempty(images{tile_id}))
      continue;
    end
    if(sift_config.is_inverted_display)
      images{tile_id} = 1-images{tile_id};
    end
    load([sift_dir, image_prefixes{tile_id}, '.sift_transform.in_plane.mat'], 'transform_p');
    transform = [transform_g([1 3 5; 2 4 6]); 0 0 1] ...
      * [transform_p([1 3 5; 2 4 6]); 0 0 1];
    transform = transform([1 2 4 5 7 8]);
    tforms{layer_id, tile_id} = maketform('affine', reshape(transform, [2 3])');
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

for layer_id = 1:size(images_t,1)
  tile_image = zeros(maxy, maxx, 3);
  t = 0;
  for tile_id = 1:size(images_t,2)
    if(isempty(images_t{layer_id, tile_id}))
      continue;
    end
    t = t+1;
    tile_image(ydata{layer_id, tile_id}(1):ydata{layer_id, tile_id}(2), ...
      xdata{layer_id, tile_id}(1):xdata{layer_id, tile_id}(2), mod(t,3)+1) = ...
      max(tile_image(ydata{layer_id, tile_id}(1):ydata{layer_id, tile_id}(2), ...
      xdata{layer_id, tile_id}(1):xdata{layer_id, tile_id}(2), mod(t,3)+1), ...
      images_t{layer_id, tile_id});
  end
  figure; imshow(tile_image);
end

return
end
