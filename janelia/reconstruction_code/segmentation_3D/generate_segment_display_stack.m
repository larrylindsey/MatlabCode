function generate_segment_display_stack(config)
% generate_segment_display_stack(config)
% Blend a color map of labelling with the grayscale image stack for
% display.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v0  03142009  init code
%

global code_dir

display_config = config.display_segmentation_3D;
if(~isfield(display_config, 'is_verbose'))
  display_config.is_verbose = true;
end
if(~isfield(display_config, 'is_verbose_figures'))
  display_config.is_verbose_figures = true;
end
if(display_config.is_verbose)
  fprintf('\nSTART: generate_segment_display_stack\n');
end

image_stack_file_name = [get_image_stack_dir(config), 'image_stack', ...
  get_image_substack_prefix(config), '.tif'];
if(display_config.is_verbose)
  fprintf('Sub-stack name:%s\n', image_stack_file_name);
end
if(exist(image_stack_file_name, 'file')~=2)
  error(['Non-existent image sub-stack ', image_stack_file_name]);
end

initial_seg_dir = [get_reconstruction_dir(config), config.segmentation_3D.dir, ...
  display_config.seg_method, '/'];
seg_file_name = [initial_seg_dir, 'seg_stack', ...
  get_image_substack_prefix(config), display_config.seg_suffix, '.raw'];
if(display_config.is_verbose)
  fprintf('Initial segmentation name:%s\n', seg_file_name);
end
if(exist(seg_file_name, 'file')~=2)
  error(['Non-existent initial segmentation ', seg_file_name]);
end

seg_mapping_dir = [get_reconstruction_dir(config), config.segmentation_3D.dir, ...
  display_config.seg_mapping_method, '/'];
seg_mapping_file_name = [seg_mapping_dir, 'seg_stack', ...
  get_image_substack_prefix(config), display_config.seg_suffix, ...
  display_config.seg_mapping_suffix, '.raw'];
if(display_config.is_verbose)
  fprintf('Initial segmentation name:%s\n', seg_mapping_file_name);
end
if(exist(seg_mapping_file_name, 'file')~=2)
  error(['Non-existent initial segmentation ', seg_mapping_file_name]);
end

if(exist(display_config.save_file_name, 'file')==2)
  warning('SEGMENT_DISPLAY_STACK:OVERWRITE', ...
    ['Overwriting a pre-existing file ', display_config.save_file_name]);
  delete(display_config.save_file_name);
end

system([code_dir, 'segmentation_3D/blend_stack_segment_b ', ...
  image_stack_file_name, ' ', seg_file_name, ' ', seg_mapping_file_name, ...
  ' ', display_config.save_file_name]);

if(display_config.is_verbose)
  fprintf('\nSTOP: superpixel_2_segment_3D_ladder\n');
end
return;
end
