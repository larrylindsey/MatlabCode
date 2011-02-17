function generate_identity_align_segment_mappings(config)
% generate_identity_align_segment_mappings(config)
% Generate identity mappings for segmentation alignment to maintain segment
% identities. Useful for combining proofread volumes.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  05052009  init. code
%

stack_config = config.stack;
align_seg_config = config.align_segmentation;
seg_config = config.segmentation_2D;

align_seg_dir = get_align_seg_dir(config);

if(~isfield(align_seg_config, 'is_verbose'))
  align_seg_config.is_verbose = true;
end
if(align_seg_config.is_verbose)
  fprintf('Aligning segment maps of multiple tiles ..\n');
  fprintf('Simply generating identity mappings.\n');
end

for case_id = stack_config.case_ids
  if(align_seg_config.is_verbose)
    fprintf('plane:%d:\n', case_id);
  end
  [image_prefixes, image_sub_dirs, is_to_be_processed] = ...
    get_image_prefixes_subdirs(config, case_id);

  if(max(is_to_be_processed)==0)
    if(align_seg_config.is_verbose)
      fprintf('All tiles in this section have is_to_be_processed=false, skipping this section\n');
    end
    continue;
  end
  
  image_set_string = get_set_string(image_prefixes);
  label_mappings = {};
  for tile_id = 1:length(image_prefixes)
    [seg_method, seg_suffix] = ...
      get_segmentation_suffix(config, image_prefixes{tile_id});
    seg_dir = [get_reconstruction_dir(config), seg_config.dir, ...
      seg_method, '/'];
    seg = load2([seg_dir, image_prefixes{tile_id}, seg_suffix, '.mat']);
    if(isfield(seg, 'label_map'))
      label_mappings{tile_id} = 0:max(seg.label_map(:)); %#ok<AGROW>
    end
    if(isfield(seg, 'locked_labels'))
      label_mappings{tile_id}(seg.locked_labels+1) = ...
        seg.locked_labels; %#ok<AGROW>
    end
  end
  
  if(align_seg_config.is_verbose)
    fprintf('Saving ... ');
  end
  save_dir = [align_seg_dir, image_sub_dirs{1}];
  check_for_dir(save_dir);
  save2([save_dir, 'sec.', num2str(case_id), '.aligned_seg_mapping', ...
    image_set_string, config.segmentation_choose.choice.seg_suffix, ...
    '.mat'], 'label_mappings');
  if(align_seg_config.is_verbose)
    fprintf('done\n');
  end
end

return
end
