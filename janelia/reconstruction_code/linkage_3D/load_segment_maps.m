function segment_2D_label_map = load_segment_maps(image_prefixes, size_images, config)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

stack_config = config.stack;
seg_config = config.segmentation_2D(end);
linkage_config = config.linkage;

segment_2D_label_map = {};
for tile = 1:length(image_prefixes)
  [seg_method, seg_suffix] = get_segmentation_suffix(config, image_prefixes{tile});
  seg_dir = [get_reconstruction_dir(config), seg_config.dir, seg_method, '/'];
  segment_2D_label_map{tile} = load2([seg_dir, image_prefixes{tile}, ...
    seg_suffix,'.mat']); %#ok<AGROW>
  if(~isfield(segment_2D_label_map{tile}, 'label_map'))
    if(isfield(config, 'segmentation_2D'))
      config_segmentation_2D = config.segmentation_2D;
      replace_flag = true;
    end
    config.segmentation_2D = config.superpixel_2_seg(1);
    [sp_method, sp_suffix] = ...
      get_superpixel_suffixes(config, image_prefixes{tile});
    if(replace_flag)
      config.segmentation_2D = config_segmentation_2D;
    end
    sp_suffix = sp_suffix{1};
    sp_dir = [get_reconstruction_dir(config), seg_config.dir, sp_method, '/'];
    superpixel_label_map = load2([sp_dir, image_prefixes{tile}, ...
      sp_suffix,'.mat']);
    superpixel_label_map.label_map(superpixel_label_map.label_map<0) = 0;
    sp_to_seg_map = [];
    sp_to_seg_map(segment_2D_label_map{tile}.superpixel_to_seg_label(:,1)+1)=...
      segment_2D_label_map{tile}.superpixel_to_seg_label(:,2); %#ok<AGROW>
    if(~isempty(sp_to_seg_map))
      segment_2D_label_map{tile}.label_map = ...
        sp_to_seg_map(1+superpixel_label_map.label_map); %#ok<AGROW>
    else
      segment_2D_label_map{tile}.label_map = []; %#ok<AGROW>
    end
  end
  segment_2D_label_map{tile}.original_label_map = segment_2D_label_map{tile}.label_map; %#ok<AGROW>
  if(isfield(linkage_config, 'modified_segment'))
    [seg_method, seg_suffix] = get_segmentation_suffix(config, image_prefixes{tile}); %#ok<ASGLU>
    seg_dir = [get_reconstruction_dir(config), seg_config.dir, ...
      linkage_config.modified_segment.method, '/'];
    temp = load2([seg_dir, image_prefixes{tile}, ...
      linkage_config.modified_segment.prefix, seg_suffix, ...
      linkage_config.modified_segment.suffix, '.mat']);
    segment_2D_label_map{tile}.label_map = temp.label_map; %#ok<AGROW>
  end
  if(~isempty(segment_2D_label_map{tile}.label_map))
    if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
      segment_2D_label_map{tile}.label_map = imresize(segment_2D_label_map{tile}.label_map, ...
        size_images{tile}, 'nearest'); %#ok<AGROW>
      segment_2D_label_map{tile}.original_label_map = ...
        imresize(segment_2D_label_map{tile}.original_label_map, ...
        size_images{tile}, 'nearest'); %#ok<AGROW>
    end
    segment_2D_label_map{tile}.label_map(segment_2D_label_map{tile}.label_map<0) = 0; %#ok<AGROW>
    segment_2D_label_map{tile}.original_label_map(segment_2D_label_map{tile}.original_label_map<0) = 0; %#ok<AGROW>
  end
end
return
end
