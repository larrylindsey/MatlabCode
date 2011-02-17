function segmentation_ncut_dir = get_segmentation_ncut_dir(config)
% segmentation_ncut_dir = get_segmentation_ncut_dir(config)
%
% Get the normalized_cuts directory from the config
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
%

seg_config = config.segmentation_2D;
ncut_config = config.normalized_cuts;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(seg_config.dir, 'dir')~=7)
  mkdir2(seg_config.dir);
end;
cd(seg_config.dir); 
if(exist(ncut_config.dir, 'dir')~=7)
  mkdir2(ncut_config.dir);
end;
cd(prev_dir);

segmentation_ncut_dir = [get_reconstruction_dir(config), seg_config.dir, ncut_config.dir];

return
end
