function segmentation_dir = get_segmentation_3D_dir(config)
% segmentation_dir = get_segmentation_dir(config)
%
% Get the stack image directory from the config
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
%

seg_config = config.segmentation_3D;
segmentation_dir = [get_reconstruction_dir(config), seg_config.dir, seg_config.method, '/'];
check_for_dir(segmentation_dir);
return
end
