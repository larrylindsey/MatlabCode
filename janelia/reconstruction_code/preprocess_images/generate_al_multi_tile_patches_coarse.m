function generate_al_multi_tile_patches_coarse(config)
% generate_al_multi_tile_patches_coarse(config)
% generate coarse al when there are multiple tiles and patches per section
% ( trakEM xml)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
% v1  09302008  modified to multiple tile per section (trakEM xml)
%

if(config.is_verbose)
  fprintf('START: generate_al_multi_tile_patches_coarse\n');
  fprintf('Constructing coarse al without roi for preview ...\n');
  fprintf('NOTE: run 810.35 to produce final al for proofreader.\n');
end

stack_config = config.stack;

load2([get_global_stitching_param_dir(config), 'canvas.mat'], ...
  'canvas_height', 'canvas_width');

align_roi_file_name = [get_global_stitching_param_dir(config), config.stack.align.roi_file_name];
load2(align_roi_file_name, 'align_roi');

fold_dir = get_fold_dir(config);
if(~isfield(config.coarse_al, 'save_sections_separately') || ...
    ~config.coarse_al.save_sections_separately)
  save_file_name = [get_region_dir(config), config.coarse_al.file_name];
  if(exist(save_file_name, 'file')==2)
    delete(save_file_name);
  end
else
  check_for_dir([get_region_dir(config), config.coarse_al.dir]);
end

for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  fprintf('plane: %d\n', case_id);
  [images, image_prefixes]  = get_image_from_stack(config, case_id);
  load2([get_global_stitching_param_dir(config), ...
    'global_stitching_parameters.c.', ...
    num2str(case_id), '.mat'], 'xdata', 'ydata', 'tforms');
  tile_image = zeros(canvas_height, canvas_width)-1;
  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))
      continue;
    end
    fprintf('tile: %d\n', tile_id);
    file_name_prefix = [fold_dir, image_prefixes{tile_id}, '.fold_mask'];
    fold_mask = load_fold_mask(file_name_prefix, config);
    for patch_id = 1:size(xdata,2) %#ok<USENS>
      if(isempty(xdata{tile_id, patch_id}))
        continue;
      end
      fprintf('patch: %d\n', patch_id);
      image = images{tile_id};
      image(fold_mask~=patch_id) = -1;
      image_t = imtransform(image, tforms{tile_id, patch_id}, ...
        'FillValue', -1); %#ok<USENS>
      temp_image = zeros(canvas_height, canvas_width)-1;
      temp_image(...
        ydata{tile_id, patch_id}(1):ydata{tile_id, patch_id}(2), ...
        xdata{tile_id, patch_id}(1):xdata{tile_id, patch_id}(2)) = ...
        image_t; %#ok<USENS>
      temp_image(tile_image>=0) = tile_image(tile_image>=0);
      tile_image = temp_image;
    end
  end
  tile_image(tile_image<0) = 0;
  al = uint8(align_roi.mask .* tile_image(align_roi.ymin:align_roi.ymax, ...
    align_roi.xmin:align_roi.xmax)*255);
  al = imresize(al, 1/config.coarse_al.scale, 'bilinear');
  if(~isfield(config.coarse_al, 'save_sections_separately') || ...
      ~config.coarse_al.save_sections_separately)
    imwrite(al, save_file_name, 'Compression', 'none', 'WriteMode', 'append');
  else
    save_file_name = sprintf([get_region_dir(config), config.coarse_al.dir, ...
      config.coarse_al.sectionwise_file_name_format], case_id);
    imwrite(al, save_file_name, 'Compression', 'none');
  end
end

if(config.is_verbose)
  fprintf('STOP: generate_al_multi_tile_patches_coarse\n');
end
return
end
