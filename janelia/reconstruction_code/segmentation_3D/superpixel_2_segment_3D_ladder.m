function superpixel_2_segment_3D_ladder(config)
% superpixel_2_segment_3D_ladder(config)
% Compute superpixels from initial watershed segmentation using ladder
%
% The initial idea for ladder was proposed by Yuriy Mishchenko, JFRC, HHMI.
% An equivalent but simpler and more efficient implementation in C for 2D
% and 3D by Shiv Vitaladevuni, JFRC, HHMI.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v0  03112009  init code
%

global code_dir

seg_config = config.segmentation_3D;
if(strcmp(seg_config.method, 'ladder')~=1)
  error('Initial segmentation method not grayscale_ladder');
end
if(~isfield(seg_config, 'is_verbose'))
  seg_config.is_verbose = true;
end
if(~isfield(seg_config, 'is_verbose_figures'))
  seg_config.is_verbose_figures = true;
end
if(seg_config.is_verbose)
  fprintf('\nSTART: superpixel_2_segment_3D_ladder\n');
end

image_stack_file_name = [get_image_stack_dir(config), 'image_stack', ...
  get_image_substack_prefix(config), '.tif'];
if(seg_config.is_verbose)
  fprintf('Sub-stack name:%s\n', image_stack_file_name);
end
if(exist(image_stack_file_name, 'file')~=2)
  error(['Non-existent image sub-stack ', image_stack_file_name]);
end

initial_seg_dir = [get_reconstruction_dir(config), config.segmentation_3D.dir, ...
  seg_config.initial_seg.method, '/'];
initial_seg_file_name = [initial_seg_dir, 'seg_stack', ...
  get_image_substack_prefix(config), seg_config.initial_seg.suffix, '.raw'];
if(seg_config.is_verbose)
  fprintf('Initial segmentation name:%s\n', initial_seg_file_name);
end
if(exist(initial_seg_file_name, 'file')~=2)
  error(['Non-existent initial segmentation ', initial_seg_file_name]);
end

seg_dir = get_segmentation_3D_dir(config);
save_suffix = [seg_config.initial_seg.suffix, '.ld.T', ...
  num2str(seg_config.f_threshold), ...
  '_L', num2str(seg_config.area_threshold)];
save_file_name = [seg_dir, 'seg_stack', get_image_substack_prefix(config), ...
  save_suffix, '.raw'];
if(exist(save_file_name, 'file')==2)
  warning('SEGMENT_3D_WATERSHED:OVERWRITE', ...
    ['Overwriting a pre-existing file ', save_file_name]);
  delete(save_file_name);
end

system([code_dir, 'segmentation_3D/superpixel_2_segment_3D_ladder_b ', ...
  image_stack_file_name, ' ', initial_seg_file_name, ' ', ...
  num2str(seg_config.f_threshold), ' ', num2str(seg_config.area_threshold), ...
  ' ', save_file_name]);

if(seg_config.is_verbose)
  fprintf('\nSTOP: superpixel_2_segment_3D_ladder\n');
end
return;
end
