function segmentation_dir = get_segmentation_dir(config)
% segmentation_dir = get_segmentation_dir(config)
%
% Get the stack image directory from the config
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
%

seg_config = config.segmentation_2D;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(seg_config.dir, 'dir')~=7)
  mkdir2(seg_config.dir);
end;
cd(seg_config.dir); 
if(exist(seg_config.method, 'dir')~=7)
  mkdir2(seg_config.method);
end;
cd(prev_dir);

segmentation_dir = [get_reconstruction_dir(config), seg_config.dir, seg_config.method, '/'];

return
end