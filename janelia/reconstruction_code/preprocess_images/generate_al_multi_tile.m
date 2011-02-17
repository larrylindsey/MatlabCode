function generate_al_multi_tile(config)
% generate_al_multi_tile(config)
% generate al when there are multiple tiles per image (trakEM xml)
%
% al{z}'s are the EM images used for the reconstruction
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
% v1  09302008  modified to multiple tile per section (trakEM xml)
%

fprintf('Constructing al ..\n');

stack_config = config.stack;
sift_dir = [get_region_dir(config), config.SIFT.dir];

align_roi_file_name = [get_region_dir(config), config.stack.align.roi_file_name];
load2(align_roi_file_name, 'align_roi');

al = {};
case_id = stack_config.case_ids(1);
fprintf('%d ', case_id);

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
    load2([sift_dir, image_prefixes{tile_id}, '.sift_transforms.mat']);
    tforms{layer_id, tile_id} = maketform('affine', reshape(transform, [2 3])'); %#ok<AGROW>
    [images_t{layer_id, tile_id}, xdata{layer_id, tile_id}, ydata{layer_id, tile_id}] = ...
      imtransform(images{tile_id}, tforms{layer_id, tile_id}, 'nearest', 'FillValues', -1); %#ok<AGROW>
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
    xdata{layer_id, tile_id} = round(xdata{layer_id, tile_id}-minx+1); %#ok<AGROW>
    ydata{layer_id, tile_id} = round(ydata{layer_id, tile_id}-miny+1); %#ok<AGROW>
  end
end
maxx = round(maxx - minx + 1);
maxy = round(maxy - miny + 1);

for layer_id = 1:size(images_t,1)
  tile_image = zeros(maxy, maxx)-1;
  for tile_id = 1:size(images_t,2)
    if(isempty(images_t{layer_id, tile_id}))
      continue;
    end
    temp_image = zeros(maxy, maxx)-1;
    temp_image(ydata{layer_id, tile_id}(1):ydata{layer_id, tile_id}(2), ...
      xdata{layer_id, tile_id}(1):xdata{layer_id, tile_id}(2)) = images_t{layer_id, tile_id};
    temp_image(tile_image>=0) = tile_image(tile_image>=0);
    tile_image = temp_image;
  end
  al_temp = uint8(align_roi.mask .* tile_image(align_roi.ymin:align_roi.ymax, ...
    align_roi.xmin:align_roi.xmax)*255);
  if(isfield(config.proofreader.roi, 'ymin'))
    al{end+1} = al_temp(config.proofreader.roi.ymin:config.proofreader.roi.ymax, ...
      config.proofreader.roi.xmin:config.proofreader.roi.xmax); %#ok<AGROW>
  else
    al{end+1} = al_temp; %#ok<AGROW>
  end
end

output_to_proofreader_stitched_grayscale(al, config);

fprintf('\ndone.\n');

