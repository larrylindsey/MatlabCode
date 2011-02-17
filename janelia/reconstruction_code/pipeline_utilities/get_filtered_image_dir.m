function filtered_image_dir = get_filtered_image_dir(config)
% filtered_image_dir = get_filtered_image_dir(config)
%
% Get the filtered image directory from the config
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  01022009  init code
%

filtered_image_config = config.precompute.filtered_image;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(filtered_image_config.dir, 'dir')~=7)
  mkdir2(filtered_image_config.dir);
end;
cd(prev_dir);

filtered_image_dir = [get_reconstruction_dir(config), filtered_image_config.dir];

return
end
