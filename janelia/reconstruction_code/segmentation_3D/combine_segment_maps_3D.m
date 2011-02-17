function combine_segment_maps_3D(config)
% combine_segment_maps_3D(config)
% Combine segment maps of sub-stacks. Useful if stack is very large.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v0  03172009  init code
%

global code_dir

stack_config = config.stack;
combine_config = config.combine_segment_3D;
if(strcmp(combine_config.adjacency_method, 'area_overlap_v0')~=1)
  error('Adjacency method not area_overlap_v0');
end
if(strcmp(combine_config.link_method, 'many_to_many')~=1)
  error('Link method not many_to_many');
end
if(~isfield(combine_config, 'is_verbose'))
  combine_config.is_verbose = true;
end
if(~isfield(combine_config, 'is_verbose_figures'))
  combine_config.is_verbose_figures = true;
end
if(combine_config.is_verbose)
  fprintf('\nSTART: combine_segment_maps_3D\n');
end

initial_seg_dir = [get_reconstruction_dir(config), config.segmentation_3D.dir, ...
  combine_config.seg_method, '/'];
seg_mapping_dir = [get_reconstruction_dir(config), config.segmentation_3D.dir, ...
  combine_config.seg_mapping_method, '/'];
combine_dir = [get_reconstruction_dir(config), config.segmentation_3D.dir, ...
  combine_config.dir];
check_for_dir(combine_dir);
save_suffix = [combine_config.seg_suffix, combine_config.seg_mapping_suffix, ...
  '.aov0_m2m_ta', num2str(combine_config.area_overlap_threshold), ...
    '_tb', num2str(combine_config.area_overlap_norm1_threshold), ...
    '_tc', num2str(combine_config.area_overlap_norm2_threshold)];
for sub_stack_id = 1:length(stack_config.sub_stack)-1
  case_ids_old = config.stack.case_ids;

  config.stack.case_ids = stack_config.sub_stack(sub_stack_id).start:...
    stack_config.sub_stack(sub_stack_id).end;
  initial_seg_file_name_1 = [initial_seg_dir, 'seg_stack', ...
    get_image_substack_prefix(config), combine_config.seg_suffix, '.raw'];
  if(combine_config.is_verbose)
    fprintf('Initial segmentation name 1:%s\n', initial_seg_file_name_1);
  end
  if(exist(initial_seg_file_name_1, 'file')~=2)
    error(['Non-existent initial segmentation 1', initial_seg_file_name_1]);
  end
  seg_mapping_file_name_1 = [seg_mapping_dir, 'seg_stack', ...
    get_image_substack_prefix(config), combine_config.seg_suffix, ...
    combine_config.seg_mapping_suffix, '.raw'];
  if(combine_config.is_verbose)
    fprintf('Initial segmentation mapping name 2:%s\n', seg_mapping_file_name_1);
  end
  if(exist(seg_mapping_file_name_1, 'file')~=2)
    error(['Non-existent initial segmentation ', seg_mapping_file_name_1]);
  end

  save_file_name_1 = [combine_dir, 'c_seg_stack', get_image_substack_prefix(config), ...
    save_suffix, '.raw'];
  if(combine_config.is_verbose)
    fprintf('Combined segmentation mapping name 1:%s\n', save_file_name_1);
  end
  if(exist(save_file_name_1, 'file')==2)
    warning('COMBINE_SEGMENT_3D:OVERWRITE', ...
      ['Overwriting a pre-existing file ', save_file_name_1]);
    delete(save_file_name_1);
  end

  config.stack.case_ids = stack_config.sub_stack(sub_stack_id+1).start:...
    stack_config.sub_stack(sub_stack_id+1).end;
  initial_seg_file_name_2 = [initial_seg_dir, 'seg_stack', ...
    get_image_substack_prefix(config), combine_config.seg_suffix, '.raw'];
  if(combine_config.is_verbose)
    fprintf('Initial segmentation name 2:%s\n', initial_seg_file_name_2);
  end
  if(exist(initial_seg_file_name_2, 'file')~=2)
    error(['Non-existent initial segmentation 2', initial_seg_file_name_2]);
  end
  seg_mapping_file_name_2 = [seg_mapping_dir, 'seg_stack', ...
    get_image_substack_prefix(config), combine_config.seg_suffix, ...
    combine_config.seg_mapping_suffix, '.raw'];
  if(combine_config.is_verbose)
    fprintf('Initial segmentation mapping name 2:%s\n', seg_mapping_file_name_2);
  end
  if(exist(seg_mapping_file_name_2, 'file')~=2)
    error(['Non-existent initial segmentation ', seg_mapping_file_name_2]);
  end

  save_file_name_2 = [combine_dir, 'c_seg_stack', get_image_substack_prefix(config), ...
    save_suffix, '.raw'];
  if(combine_config.is_verbose)
    fprintf('Combined segmentation mapping name 2:%s\n', save_file_name_2);
  end
  if(exist(save_file_name_2, 'file')==2)
    warning('COMBINE_SEGMENT_3D:OVERWRITE', ...
      ['Overwriting a pre-existing file ', save_file_name_2]);
    delete(save_file_name_2);
  end

  n_section_overlap = stack_config.sub_stack(sub_stack_id).end - ...
    stack_config.sub_stack(sub_stack_id+1).start + 1;
  if(n_section_overlap<=0)
    error('Number of sections overlap is less than 0. Check combine_segment_3D.sub_stack');
  end
  if(combine_config.is_verbose)
    fprintf('Number of overlapping sections :%d\n', n_section_overlap);
  end
  
  system([code_dir, 'segmentation_3D/combine_segment_maps_3D_aov0_m2m_b', ...
    ' ', initial_seg_file_name_1, ' ', seg_mapping_file_name_1, ...
    ' ', initial_seg_file_name_2, ' ', seg_mapping_file_name_2, ...
    ' ', num2str(n_section_overlap), ...
    ' ', num2str(combine_config.area_overlap_threshold), ...
    ' ', num2str(combine_config.area_overlap_norm1_threshold), ...
    ' ', num2str(combine_config.area_overlap_norm2_threshold), ...
    ' ', save_file_name_1, ' ', save_file_name_2]);

  config.stack.case_ids = case_ids_old ;
end

if(combine_config.is_verbose)
  fprintf('\nSTOP: combine_segment_maps_3D\n');
end
return;
end
