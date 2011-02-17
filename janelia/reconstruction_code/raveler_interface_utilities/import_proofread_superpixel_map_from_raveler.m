function import_proofread_superpixel_map_from_raveler(config)
% import_proofread_superpixel_map_from_raveler(config)
% Imports proofread superpixel maps from files produced by Raveler,
% typically saved after using the matlab gui.m
% Copies onto the image plane the superpixel map from that shown to the
% user during proofreading (possibly stitched from multiple tiles).
% Label '0' is written on all undefined pixels.
%
% Shiv N. Vitaladevuni
% v0  04062009  init. code
%

fprintf('START: import_proofread_superpixel_map_from_raveler\n');
if(strcmp(config.proofreader.method, 'Raveler')~=1)
  error('Proofreader method is not Raveler');
end

fprintf('Importing proofread data from proofreader ..\n');

stack_config = config.stack;
import_config = config.proofreader.import;
raveler_config = import_config.Raveler;

switch(config.align.global_align.method)
  case 'SIFT'
    tform_dir = [raveler_config.to_be_proofread_dir, config.SIFT.dir];
  case 'deformable_mesh'
    tform_dir = [raveler_config.to_be_proofread_dir, config.deformable_mesh.dir];
end

load2([raveler_config.to_be_proofread_dir, 'global_stitching_parameters/', ...
  'canvas.mat'], 'canvas_height', 'canvas_width');

align_roi_file_name = [raveler_config.to_be_proofread_dir, 'global_stitching_parameters/', ...
  config.stack.align.roi_file_name];
load2(align_roi_file_name, 'align_roi');

proofreader_roi_file_name = ...
  [raveler_config.to_be_proofread_dir, 'global_stitching_parameters/', ...
  config.proofreader.roi.file_name];
try
  load2(proofreader_roi_file_name, 'proofreader_roi');
catch %#ok<CTCH>
  proofreader_roi = struct([]);
end

save_dir = [get_reconstruction_dir(config), ...
  config.superpixel(1).dir, import_config.dir];
check_for_dir(save_dir);
save_suffix = ['.prfrd_csp', '.', import_config.proofread_name];
fold_dir = [raveler_config.to_be_proofread_dir, 'fold_masks/'];
for case_id = stack_config.case_ids
  layer_id = find(config.region.case_ids==case_id);
  fprintf('plane: %d, layer_id: %d\n', case_id, layer_id);
  
  [image_prefixes, image_sub_dirs] = ...
    get_image_prefixes_subdirs(config, case_id);
  
  load2([raveler_config.to_be_proofread_dir, 'global_stitching_parameters/' ...
    'global_stitching_parameters.c.' num2str(case_id), '.mat'], ...
    'xdata', 'ydata', 'tforms');
  
  if(isfield(raveler_config, 'is_numbered_by_depth') && ...
      ~raveler_config.is_numbered_by_depth)
    file_name_import_superpixel_map = sprintf([...
      raveler_config.proofread_data_dir, raveler_config.superpixel_dir, ...
      raveler_config.superpixel_prefix, raveler_config.superpixel_suffix], ...
      layer_id);
  else
    file_name_import_superpixel_map = sprintf([...
      raveler_config.proofread_data_dir, raveler_config.superpixel_dir, ...
      raveler_config.superpixel_prefix, raveler_config.superpixel_suffix], ...
      case_id);
  end
  stitched_superpixel_map = imread(file_name_import_superpixel_map);
  
  stitched_label_map_1 = zeros(align_roi.ymax-align_roi.ymin+1, ...
    align_roi.xmax-align_roi.xmin+1);
  if(isfield(proofreader_roi, 'ymin'))
    stitched_label_map_1(proofreader_roi.ymin:proofreader_roi.ymax, ...
      proofreader_roi.xmin:proofreader_roi.xmax) = stitched_superpixel_map;
  else
    stitched_label_map_1 = stitched_superpixel_map;
  end
  stitched_label_map = zeros(canvas_height, canvas_width);
  stitched_label_map(align_roi.ymin:align_roi.ymax, align_roi.xmin:align_roi.xmax) = ...
    stitched_label_map_1;
  stitched_label_map = uint32(stitched_label_map);
  
  fold_mask = cell(1, length(image_prefixes));
  patch_id_offset = zeros(length(image_prefixes), 1);
  transforms_affine = [];
  fprintf('Generating patch_id maps ...\n');
  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))
      continue;
    end
    fprintf('tile: %d\n%s\n', tile_id, image_prefixes{tile_id});
    file_name_prefix = [fold_dir, image_prefixes{tile_id}, '.fold_mask'];
    fold_mask{tile_id} = ...
      uint32(load_fold_mask(file_name_prefix, config));
    fold_mask{tile_id}(fold_mask{tile_id}>0) = ...
      fold_mask{tile_id}(fold_mask{tile_id}>0) + patch_id_offset(tile_id);
    
%     load2([tform_dir, image_prefixes{tile_id}, ...
%       '.patchwise_affine_transforms.mat'], 'transform');
    fin = fopen([tform_dir, image_prefixes{tile_id}, ...
      '.global_patchwise_tform.txt'], 'rt');
    transform.transform = fscanf(fin, '%g', [6 inf]);
    fclose(fin);
    transform.transform = transform.transform([1 4 2 5 3 6], :);
    transforms_affine = [transforms_affine, transform.transform]; %#ok<AGROW>
    
    patch_id_offset(tile_id+1) = max(fold_mask{tile_id}(:));
  end
  fprintf('done.\n');
  display(transforms_affine)
  
  % render the fold masks in global coordinates (overlapping).
  fprintf('Rendering fold masks in global coordinates ...\n');
  section_image = zeros(canvas_height, canvas_width)-1;
  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))
      continue;
    end
    fprintf('tile: %d\n%s\n', tile_id, image_prefixes{tile_id});
    m_p_id = 0;
    for p_id = 1:size(xdata,2)
      if(isempty(xdata{tile_id, p_id}))
        continue;
      end
      m_p_id = max(m_p_id, p_id);
    end
    m_p_id_1 = max(fold_mask{tile_id}(:));
    if(m_p_id_1>0 && m_p_id_1 ~= m_p_id+patch_id_offset(tile_id))
      fprintf('maximum patch id in fold mask: %d\n', m_p_id_1);
      fprintf('maximum patch id in loaded transform: %d\n', m_p_id+patch_id_offset(tile_id));
      error('Number of patches in saved transforms doesn''t match with number of patches in fold mask.');
    end
    for patch_id = 1:size(xdata,2) %#ok<*USENS>
      if(isempty(xdata{tile_id, patch_id}))
        continue;
      end
      fprintf('patch_id: %d\n', patch_id);
      image = double(fold_mask{tile_id});
      image(image~=patch_id+patch_id_offset(tile_id)) = -1;
      [image_t, xdata_1, ydata_1] = ...
        imtransform(image, tforms{tile_id, patch_id}, ...
        'FillValue', -1); %#ok<USENS>
      temp_image = zeros(canvas_height, canvas_width)-1;
      temp_image(...
        ydata{tile_id, patch_id}(1):ydata{tile_id, patch_id}(2), ...
        xdata{tile_id, patch_id}(1):xdata{tile_id, patch_id}(2)) = ...
        image_t; %#ok<USENS>
      temp_image(section_image>=0) = section_image(section_image>=0);
      section_image = temp_image;
      transforms_affine(5,patch_id+patch_id_offset(tile_id)) = ...
        transforms_affine(5,patch_id+patch_id_offset(tile_id)) - ...
        xdata_1(1) + xdata{tile_id, patch_id}(1); %#ok<AGROW>
      transforms_affine(6,patch_id+patch_id_offset(tile_id)) = ...
        transforms_affine(6,patch_id+patch_id_offset(tile_id)) - ...
        ydata_1(1) + ydata{tile_id, patch_id}(1); %#ok<AGROW>
    end
  end
  display(transforms_affine)
  section_image(section_image<0) = 0;
  section_image(1:align_roi.ymin-1, :) = 0;
  section_image(align_roi.ymax+1:end, :) = 0;
  section_image(:, 1:align_roi.xmin-1) = 0;
  section_image(:, align_roi.xmax+1:end) = 0;
  section_image(align_roi.ymin:align_roi.ymax, align_roi.xmin:align_roi.xmax) = ...
    align_roi.mask .* ...
    section_image(align_roi.ymin:align_roi.ymax, align_roi.xmin:align_roi.xmax);
  section_image = uint32(section_image);
  
  if(isfield(proofreader_roi, 'ymin'))
    xmin = align_roi.xmin + proofreader_roi.xmin - 1;
    ymin = align_roi.ymin + proofreader_roi.ymin - 1;
    xmax = align_roi.xmin + proofreader_roi.xmax - 1;
    ymax = align_roi.ymin + proofreader_roi.ymax - 1;
    section_image(1:ymin-1, :) = 0;
    section_image(ymax+1:end, :) = 0;
    section_image(:, 1:xmin-1) = 0;
    section_image(:, xmax+1:end) = 0;
  end
  fprintf('done.\n');
  
  fprintf('Copying superpixel maps from proofread onto image planes ...\n');
  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))  
      continue;
    end
    fprintf('tile: %d\n%s\n', tile_id, image_prefixes{tile_id});
    tile_image = zeros(canvas_height, canvas_width)-1;
    fprintf('Rendering patches in global coordinates ...\n');
    for patch_id = 1:size(xdata,2)
      if(isempty(xdata{tile_id, patch_id}))
        continue;
      end
      fprintf('patch_id: %d\n', patch_id);
      image = double(fold_mask{tile_id});
      image(image~=patch_id+patch_id_offset(tile_id)) = -1;
      image_t = imtransform(image, tforms{tile_id, patch_id}, ...
        'FillValue', -1);
      temp_image = zeros(canvas_height, canvas_width)-1;
      temp_image(...
        ydata{tile_id, patch_id}(1):ydata{tile_id, patch_id}(2), ...
        xdata{tile_id, patch_id}(1):xdata{tile_id, patch_id}(2)) = ...
        image_t;
      temp_image(tile_image>=0) = tile_image(tile_image>=0);
      tile_image = temp_image;
    end
    fprintf('done.\n');
    
    fprintf('Applying ROI ... ');
    tile_image(tile_image<0) = 0;
    tile_image(1:align_roi.ymin-1, :) = 0;
    tile_image(align_roi.ymax+1:end, :) = 0;
    tile_image(:, 1:align_roi.xmin-1) = 0;
    tile_image(:, align_roi.xmax+1:end) = 0;
    tile_image(align_roi.ymin:align_roi.ymax, align_roi.xmin:align_roi.xmax) = ...
      align_roi.mask .* ...
      tile_image(align_roi.ymin:align_roi.ymax, align_roi.xmin:align_roi.xmax);
    tile_image = uint32(tile_image);
    
    if(isfield(proofreader_roi, 'ymin'))
      xmin = align_roi.xmin + proofreader_roi.xmin - 1;
      ymin = align_roi.ymin + proofreader_roi.ymin - 1;
      xmax = align_roi.xmin + proofreader_roi.xmax - 1;
      ymax = align_roi.ymin + proofreader_roi.ymax - 1;
      tile_image(1:ymin-1, :) = 0;
      tile_image(ymax+1:end, :) = 0;
      tile_image(:, 1:xmin-1) = 0;
      tile_image(:, xmax+1:end) = 0;
    end
    fprintf('done.\n');
    
    % copy segmentation labels to image plane from the global coordinates
    % based on the stitched fold masks
    fprintf('Copying superpixel maps from proofread onto image planes ...\n');
    imported_superpixel_map = get_label_piecewise_affine_transformation_mask(...
      double(stitched_label_map), uint32(fold_mask{tile_id}), ...
      double(transforms_affine), -1, uint32(tile_image));
    fprintf('done.\n');
    
    tile_overlap_map = get_label_piecewise_affine_transformation_mask(...
      double(section_image), uint32(fold_mask{tile_id}), ...
      double(transforms_affine), -1, uint32(tile_image));
    
    tile_no_overlap_mask = ismember(tile_overlap_map, ...
      [0, patch_id_offset(tile_id)+1:patch_id_offset(tile_id+1)]);
    label_map = zeros(size(imported_superpixel_map))-1;
    label_map(tile_no_overlap_mask) = ...
      imported_superpixel_map(tile_no_overlap_mask);
    fprintf('done.\n');
    
    fprintf('Saving ... ');
    check_for_dir([save_dir, image_sub_dirs{tile_id}]);
    file_name = [save_dir, image_prefixes{tile_id}, ...
      save_suffix, '.mat'];
    locked_labels = unique(label_map(label_map>0)); %#ok<NASGU>
    save2(file_name, 'label_map', 'locked_labels');
    fprintf('done.\n');
  end
end

fprintf('STOP: import_proofread_superpixel_map_from_raveler\n');
return
end
