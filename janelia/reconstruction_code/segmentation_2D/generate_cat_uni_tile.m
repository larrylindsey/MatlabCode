function generate_cat_uni_tile(config)
% generate_cat_uni_tile(config)
% generate cat - one tile per section (non trakEM xml)
%
% cat{z}'s are the superpixel maps being used in the reconstruction
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
%

fprintf('Constructing cat ..\n');

stack_config = config.stack;
superpixel_config = config.superpixel(end);

image_dir = get_stack_dir(config);

% stack is not prealigned then compute an ROI for the aligned stack if
% needed
if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned ...
    && isempty(stack_config.align.roi))
  align_roi_file_name = [image_dir, stack_config.align.roi_file_name];
  if(exist(align_roi_file_name, 'file')~=2)
    align_roi = get_roi(config);
    save2(align_roi_file_name, 'align_roi');
  else
    load2(align_roi_file_name, 'align_roi');
  end;
end

case_id = stack_config.case_ids(1);
fprintf('%d ', case_id);
cat = [];

[images, image_prefixes] = get_image_from_stack(config, case_id);
image = images{1};
image_prefix = image_prefixes{1};
[superpixel_method, superpixel_suffix] = ...
  get_superpixel_suffixes_private(config, image_prefix);
superpixel_dir = [get_reconstruction_dir(config), superpixel_config.dir, ...
  superpixel_method, '/'];
superpixel_seg = load2([superpixel_dir, image_prefix, superpixel_suffix, '.mat']); 
if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
  superpixel_seg.label_map = imresize(superpixel_seg.label_map, ...
    size(image), 'nearest');
end

if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned)
  align_config = stack_config.align;
  align_xg_file_name = [image_dir, align_config.xg_file_name];
  xg = read_xf(align_xg_file_name);
  save2([image_dir, align_config.xg_save_name], 'xg');
  
  margin.left = 1-align_config.margin;
  margin.right = size(superpixel_seg.label_map,2)+align_config.margin;
  margin.top = 1-align_config.margin;
  margin.bottom = size(superpixel_seg.label_map,1)+align_config.margin;

  align_tform = xf_2_tform(xg(1,:), size(superpixel_seg.label_map, 1), ...
    size(superpixel_seg.label_map, 2));
  superpixel_seg.label_map = apply_tform(superpixel_seg.label_map, ...
    align_tform, margin, 'nearest');
  superpixel_seg.label_map = superpixel_seg.label_map(align_roi.ymin:align_roi.ymax, ...
    align_roi.xmin:align_roi.xmax);
end;

proofread_config = config.proofreader;
if(isfield(proofread_config.roi, 'ymin'))
  cat{end+1} = superpixel_seg.label_map(...
    proofread_config.roi.ymin:proofread_config.roi.ymax, ...
    proofread_config.roi.xmin:proofread_config.roi.xmax);
else
  cat{end+1} = superpixel_seg.label_map;
end


for i = 2:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  fprintf('%d ', case_id);
  if(mod(i,20)==0)
    fprintf('\n');
  end;
  
  [images, image_prefixes] = get_image_from_stack(config, case_id);
  image = images{1};
  image_prefix = image_prefixes{1};
  [superpixel_method, superpixel_suffix] = ...
    get_superpixel_suffixes_private(config, image_prefix);
  superpixel_dir = [get_reconstruction_dir(config), superpixel_config.dir, ...
    superpixel_method, '/'];
  superpixel_seg = load2([superpixel_dir, image_prefix, superpixel_suffix, '.mat']);
  if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
    superpixel_seg.label_map = imresize(superpixel_seg.label_map, ...
      size(image), 'nearest');
  end

  if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned)
    margin.left = 1-align_config.margin;
    margin.right = size(superpixel_seg.label_map,2)+align_config.margin;
    margin.top = 1-align_config.margin;
    margin.bottom = size(superpixel_seg.label_map,1)+align_config.margin;

    align_tform = xf_2_tform(xg(i,:), size(superpixel_seg.label_map, 1), ...
      size(superpixel_seg.label_map, 2));
    superpixel_seg.label_map = apply_tform(superpixel_seg.label_map, ...
      align_tform, margin, 'nearest');
    superpixel_seg.label_map = superpixel_seg.label_map(align_roi.ymin:align_roi.ymax, ...
      align_roi.xmin:align_roi.xmax);
  end;
  
  if(isfield(proofread_config.roi, 'ymin'))
    cat{end+1} = superpixel_seg.label_map(...
      proofread_config.roi.ymin:proofread_config.roi.ymax, ...
      proofread_config.roi.xmin:proofread_config.roi.xmax);
  else
    cat{end+1} = superpixel_seg.label_map;
  end
end;

output_to_proofreader_stitched_superpixel_maps(cat, config);

fprintf('\ndone.\n');
return

  function [superpixel_method, superpixel_suffix] = ...
    get_superpixel_suffixes_private(config, image_prefix)
    replace_flag = false;
    if(isfield(config, 'segmentation_2D'))
      config_segmentation_2D = config.segmentation_2D;
    end
    config.segmentation_2D = config.superpixel_2_seg(1);
    [superpixel_method, superpixel_suffixes] = ...
      get_superpixel_suffixes(config, image_prefix);
    superpixel_suffix = superpixel_suffixes{1};
    if(replace_flag)
      config.segmentation_2D = config_segmentation_2D;
    end
    return
  end
end
