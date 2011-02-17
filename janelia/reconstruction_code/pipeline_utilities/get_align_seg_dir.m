function align_seg_dir = get_align_seg_dir(config)
% align_seg_dir = get_align_seg_dir(config)
%
% Get the segment map alignment saving directory from the config
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
%

align_config = config.align_segmentation;
seg_config = config.segmentation_2D;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(seg_config.dir, 'dir')~=7)
  mkdir2(seg_config.dir);
end;
cd(seg_config.dir); 
if(exist(align_config.dir, 'dir')~=7)
  mkdir2(align_config.dir);
end;
cd(prev_dir);

align_seg_dir = [get_reconstruction_dir(config), seg_config.dir, align_config.dir];

return
end
