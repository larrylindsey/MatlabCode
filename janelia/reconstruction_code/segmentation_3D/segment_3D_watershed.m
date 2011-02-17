function segment_3D_watershed(config)
% segment_3D_watershed(config)
% Perform 3D segmentation using watershed
%
% Uses watershed implementation by Gene Myers, JFRC, HHMI.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v0  03112009  init code
%

global code_dir

seg_config = config.segmentation_3D;
if(strcmp(seg_config.method, 'watershed')~=1)
  error('Initial segmentation method not watershed');
end
if(~isfield(seg_config, 'is_verbose'))
  seg_config.is_verbose = true;
end
if(~isfield(seg_config, 'is_verbose_figures'))
  seg_config.is_verbose_figures = true;
end
if(seg_config.is_verbose)
  fprintf('\nSTART: segment_3D_watershed\n');
end

image_stack_file_name = [get_image_stack_dir(config), 'image_stack', ...
  get_image_substack_prefix(config), '.tif'];
if(seg_config.is_verbose)
  fprintf('Sub-stack name:%s\n', image_stack_file_name);
end
if(exist(image_stack_file_name, 'file')~=2)
  error(['Non-existent image sub-stack file ', image_stack_file_name]);
end

seg_dir = get_segmentation_3D_dir(config);
save_suffix = ['.ws.T', num2str(seg_config.f_threshold)];
save_file_name = [seg_dir, 'seg_stack', get_image_substack_prefix(config), ...
  save_suffix, '.raw'];
if(exist(save_file_name, 'file')==2)
  warning('SEGMENT_3D_WATERSHED:OVERWRITE', ...
    ['Overwriting a pre-existing file ', save_file_name]);
  delete(save_file_name);
end

system([code_dir, 'segmentation_3D/watershed_3D_myers ', ...
  image_stack_file_name, ' ', num2str(seg_config.f_threshold), ' ', save_file_name]);

if(seg_config.is_verbose)
  fprintf('\nSTOP: segment_3D_watershed\n');
end
return;
end
