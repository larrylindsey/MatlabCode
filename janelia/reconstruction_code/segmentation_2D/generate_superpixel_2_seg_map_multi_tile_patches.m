function generate_superpixel_2_seg_map_multi_tile_patches(config)
% generate_superpixel_2_seg_map_multi_tile_patches(config)
% generate superpixel_2_seg_map - multiple tiles and patches per
% section.
%
% superpixel_2_seg_map{z}'s are the superpixel to segmentation mappings.
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('START: generate_superpixel_2_seg_map_multi_tile_patches\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stack_config = config.stack;
seg_config = config.superpixel_2_seg(end);

align_seg_dir = get_align_seg_dir(config);

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
  fprintf('Merging the superpixel_2_seg_maps\n');
  
  superpixel_2_seg_map = [];

  [image_prefixes, image_sub_dirs] = ...
    get_image_prefixes_subdirs(config, case_id);

  image_set_string = get_set_string(image_prefixes);
  save_dir = [align_seg_dir, image_sub_dirs{1}];
  label_mappings = cell(1, length(image_prefixes));
  load2([save_dir, 'sec.', num2str(case_id), ...
    '.aligned_seg_mapping', image_set_string, ...
    config.segmentation_choose.choice.seg_suffix, '.mat'], 'label_mappings');
  
  % read in superpixel remapping file
  superpixel_mappings = get_segment_label_remapping([get_region_dir(config), ...
    sprintf(export_config.superpixel_remapping_file, case_id)]);
  
  superpixel_ids_sets = [];
  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))
      continue;
    end
    fprintf('tile: %d\nimage_prefix: %s\n', tile_id, image_prefixes{tile_id});

    % for the superpixel_to_seg_label ..
    % load it
    [seg_method, seg_suffix] = ...
      get_segmentation_suffix(config, image_prefixes{tile_id});
    seg_dir = [get_reconstruction_dir(config), seg_config.dir, seg_method, '/'];
    seg_label_map = load2([seg_dir, image_prefixes{tile_id}, ...
      seg_suffix,'.mat']);
    
    % find the relevant superpixel mapping
    tt = -1;
    for t = 1:length(superpixel_mappings)
      if(~isempty(strfind(superpixel_mappings(t).segment_file_name, image_prefixes{tile_id})))
        tt = t;
        break;
      end
    end
    if(tt==-1)
      error('Could not find label mapping for this tile');
    end
    
    if(~isempty(seg_label_map.superpixel_to_seg_label))
      % add superpixel label offset
      seg_label_map.superpixel_to_seg_label(:,1) = ...
        apply_mapping(seg_label_map.superpixel_to_seg_label(:,1), ...
        superpixel_mappings(tt).label_mapping);
      % include mapping of segment ids from the stitching
      label_mappings{tile_id}(1) = 0; % 0 label always gets mapped to 0
      seg_label_map.superpixel_to_seg_label(:,2) = ...
        label_mappings{tile_id}(seg_label_map.superpixel_to_seg_label(:,2)+1);
      
      % remove invalid superpixel labels
      seg_label_map.superpixel_to_seg_label = ...
        seg_label_map.superpixel_to_seg_label(...
        seg_label_map.superpixel_to_seg_label(:,1)>0, :);
      
      % update the layer's superpixel_2_seg_map
      superpixel_2_seg_map(seg_label_map.superpixel_to_seg_label(:,1)+1) = ...
        seg_label_map.superpixel_to_seg_label(:,2); %#ok<AGROW>
      
      superpixel_ids_sets = ...
        [superpixel_ids_sets, seg_label_map.superpixel_to_seg_label(:,1)']; %#ok<AGROW>
    end
  end
  
  if(~isempty(superpixel_2_seg_map))
    % due to upscaling it is possible that some superpixels which were
    % covered during stitching may become exposed. So set any negative
    % superpixel_2_seg_map values to 0
    superpixel_2_seg_map(superpixel_2_seg_map<0) = 0; %#ok<AGROW>
  end
  
  superpixel_2_seg_map_t = superpixel_2_seg_map; %#ok<NASGU>
  save2([get_region_dir(config), 'superpixel_2_seg_map_t.', num2str(case_id), ...
    config.segmentation_choose.choice.seg_suffix, '.mat'], ...
    'superpixel_2_seg_map_t', 'superpixel_ids_sets');
end
fprintf('STOP: generate_superpixel_2_seg_map_multi_tile_patches\n');

return
end
