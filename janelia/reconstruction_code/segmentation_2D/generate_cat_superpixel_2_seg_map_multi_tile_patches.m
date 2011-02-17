function generate_cat_superpixel_2_seg_map_multi_tile_patches(config)
% generate_cat_superpixel_2_seg_map_multi_tile_patches(config)
% generate cat and superpixel_2_seg_map - multiple tiles and patches per
% section.
%
% cat{z}'s are the superpixel maps being used in the reconstruction
% superpixel_2_seg_map{z}'s are the superpixel to segmentation mappings.
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
% v1  09302008  modified for multiple tiles per section (trakEM xml).
%                 Including generation of superpixel_2_seg_map.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Constructing cat and superpixel_2_seg_map ..\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stack_config = config.stack;
seg_config = config.superpixel_2_seg(end);
superpixel_config = config.superpixel(end);

align_seg_dir = get_align_seg_dir(config);

load2([get_global_stitching_param_dir(config), 'canvas.mat'], ...
  'canvas_height', 'canvas_width');

align_roi_file_name = [get_global_stitching_param_dir(config), ...
  stack_config.align.roi_file_name];
load2(align_roi_file_name, 'align_roi');

fold_dir = get_fold_dir(config);
images = get_image_from_stack(config, stack_config.case_ids(1));
size_image = size(images{1});

proofread_config = config.proofreader;
if(~isfield(proofread_config, 'export'))
  proofread_config.export = struct([]);
end
export_config = proofread_config.export;
if(~isfield(export_config, 'is_verbose_figures'))
  export_config.is_verbose_figures = false;
end
  
for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  fprintf('plane: %d\n', case_id);
  fprintf('Stitching the superpixel maps and merging the superpixel_2_seg_maps.\n');
  
  superpixel_2_seg_map = [];

  [image_prefixes, image_sub_dirs] = ...
    get_image_prefixes_subdirs(config, case_id);
  load2([get_global_stitching_param_dir(config), ...
    'global_stitching_parameters.c.', ...
    num2str(case_id), '.mat'], 'xdata', 'ydata', 'tforms');

  image_set_string = get_set_string(image_prefixes);
  save_dir = [align_seg_dir, image_sub_dirs{1}];
  label_mappings = cell(1, length(image_prefixes));
  load2([save_dir, 'sec.', num2str(case_id), ...
    '.aligned_seg_mapping', image_set_string, ...
    config.segmentation_choose.choice.seg_suffix, '.mat'], 'label_mappings');
  
  superpixel_id_offset = 0;
  cat_l = zeros(canvas_height, canvas_width) - 1;
  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))
      continue;
    end
    fprintf('tile: %d\n', tile_id);
    % load and transform the superpixel segmentation map
    replace_flag = false;
    if(isfield(config, 'segmentation_2D'))
      config_segmentation_2D = config.segmentation_2D;
      replace_flag = true;
    end
    config.segmentation_2D = config.superpixel_2_seg(1);
    [superpixel_method, superpixel_suffixes] = ...
      get_superpixel_suffixes(config, image_prefixes{tile_id});
    if(replace_flag)
      config.segmentation_2D = config_segmentation_2D;
    end
    superpixel_dir = [get_reconstruction_dir(config), superpixel_config.dir, ...
      superpixel_method, '/'];
    superpixel_seg = load2([superpixel_dir, image_prefixes{tile_id}, ...
      superpixel_suffixes{1}, '.mat']);
    if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
      superpixel_seg.label_map = imresize(superpixel_seg.label_map, ...
        size_image, 'nearest');
    end
    
    file_name_prefix = ...
        [fold_dir, image_prefixes{tile_id}, '.fold_mask'];
    fold_mask = load_fold_mask(file_name_prefix, config);
    temp_superpixel_map = zeros(canvas_height, canvas_width) - 1;
    for patch_id = 1:size(xdata,2) %#ok<USENS>
      if(isempty(xdata{tile_id, patch_id}))
        continue;
      end
      fprintf('patch_id: %d\n', patch_id);
      xd = xdata{tile_id, patch_id};
      yd = ydata{tile_id, patch_id}; %#ok<USENS>
      label_map_temp = superpixel_seg.label_map;
      label_map_temp(fold_mask~=patch_id) = -1;
      seg_temp = imtransform(label_map_temp, tforms{tile_id, patch_id}, ...
        'nearest', 'FillValues', -1); %#ok<USENS>
      temp_superpixel_map(yd(1):yd(2), xd(1):xd(2)) = ...
        max(temp_superpixel_map(yd(1):yd(2), xd(1):xd(2)), seg_temp);
      if(export_config.is_verbose_figures)
        figure;
        imshow(imresize(label_map_temp, 1/16, 'nearest')>0);
        title(['label_map_temp case_id: ', num2str(case_id), ...
          ' tile_id: ', num2str(tile_id), 'patch_id: ', num2str(patch_id)]);
        figure;
        imshow(imresize(temp_superpixel_map, 1/16, 'nearest'));
        title(['temp_superpixel_map case_id: ', num2str(case_id), ...
          ' tile_id: ', num2str(tile_id), 'patch_id: ', num2str(patch_id)]);
      end
    end
    % add an offset to the superpixel labels to ensure uniqueness among
    % tiles
    temp_superpixel_map(temp_superpixel_map>0) = ...
      temp_superpixel_map(temp_superpixel_map>0) + superpixel_id_offset;
    temp_superpixel_map(cat_l>=0) = cat_l(cat_l>=0);
    cat_l = temp_superpixel_map;
    if(export_config.is_verbose_figures)
      figure;
      imshow(imresize(cat_l, 1/16, 'nearest'));
      title(['cat_l case_id: ', num2str(case_id), ' tile_id: ', num2str(tile_id)]);
    end
    % for the superpixel_to_seg_label ..
    % load it
    [seg_method, seg_suffix] = ...
      get_segmentation_suffix(config, image_prefixes{tile_id});
    seg_dir = [get_reconstruction_dir(config), seg_config.dir, seg_method, '/'];
    seg_label_map = load2([seg_dir, image_prefixes{tile_id}, ...
      seg_suffix,'.mat']);
    
    if(~isempty(seg_label_map.superpixel_to_seg_label))
      % add superpixel label offset
      seg_label_map.superpixel_to_seg_label(:,1) = ...
        seg_label_map.superpixel_to_seg_label(:,1) + superpixel_id_offset;
      % include mapping of segment ids from the stitching
      label_mappings{tile_id}(1) = 0; % 0 label always gets mapped to 0
      seg_label_map.superpixel_to_seg_label(:,2) = ...
        label_mappings{tile_id}(seg_label_map.superpixel_to_seg_label(:,2)+1);
      % update the layer's superpixel_2_seg_map
      superpixel_2_seg_map(seg_label_map.superpixel_to_seg_label(:,1)+1) = ...
        seg_label_map.superpixel_to_seg_label(:,2); %#ok<AGROW>
    end
    
    % maximum superpixel id seen so far
    superpixel_id_offset = max([0, max(cat_l(:))]);
  end
  
  fprintf('Applying ROI for cat and making superpixel ids unique within sections\n');
  cat_l(cat_l<0) = 0;
  cat_l = align_roi.mask .* cat_l(...
    align_roi.ymin:align_roi.ymax, align_roi.xmin:align_roi.xmax);
  
  if(isfield(proofread_config.roi, 'ymin'))
    cat_l = cat_l(...
      proofread_config.roi.ymin:proofread_config.roi.ymax, ...
      proofread_config.roi.xmin:proofread_config.roi.xmax);
  end
  
  % Remove superpixel ids masked out by roi_mask. Relabel to remove
  % redundant labels.
  superpixel_ids_unique = unique(cat_l(:));
  superpixel_id_map = [];
  superpixel_id_map(1) = 0;
  superpixel_id_map(superpixel_ids_unique+1) = ...
    0:(length(superpixel_ids_unique)-1); %#ok<AGROW>
  cat_l = superpixel_id_map(cat_l + 1);
  
  % Some very small superpixels may appear in the global alignment that were
  % invisible in the within-section stitching. Set these to 0.
  new_superpixel_ids = superpixel_ids_unique(superpixel_ids_unique >= ...
    length(superpixel_2_seg_map));
  superpixel_2_seg_map(new_superpixel_ids+1) = 0; %#ok<AGROW>
  if(~isempty(superpixel_2_seg_map))
    superpixel_2_seg_map = ...
      superpixel_2_seg_map(superpixel_ids_unique+1);

    % due to upscaling it is possible that some superpixels which were
    % covered during stitching may become exposed. So set any negative
    % superpixel_2_seg_map values to 0
    superpixel_2_seg_map(superpixel_2_seg_map<0) = 0;
    
  end
  
  if(~isempty(superpixel_id_map))
    superpixel_ids_sets = superpixel_id_map(superpixel_ids_unique+1); %#ok<NASGU>
  else
    superpixel_ids_sets = []; %#ok<NASGU>
  end

  fprintf('Dumping cat and superpixel_2_seg_map\n');
  output_to_raveler_stitched_superpixel_map(cat_l, case_id, config);

  superpixel_2_seg_map_t = superpixel_2_seg_map; %#ok<NASGU>
  save2([get_region_dir(config), 'superpixel_2_seg_map_t.', num2str(case_id), ...
    config.segmentation_choose.choice.seg_suffix, '.mat'], ...
    'superpixel_2_seg_map_t', 'superpixel_ids_sets');
end
fprintf('done\n');

return
end
