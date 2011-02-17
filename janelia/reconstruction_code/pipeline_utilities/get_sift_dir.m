function sift_dir = get_sift_dir(config)
% sift_dir = get_sift_dir(config)
%
% Get the sift directory from the config
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
%

sift_config = config.SIFT;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(sift_config.dir, 'dir')~=7)
  mkdir2(sift_config.dir);
end;
cd(prev_dir);

sift_dir = [get_reconstruction_dir(config), sift_config.dir];

return
end
