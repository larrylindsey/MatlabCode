function linkage_dir = get_linkage_dir(config)
% linkage_dir = get_linkage_dir(config)
%
% Get the stack image directory from the config
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
%

linkage_config = config.linkage;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(linkage_config.dir, 'dir')~=7)
  mkdir2(linkage_config.dir);
end;
cd(prev_dir);

linkage_dir = [get_reconstruction_dir(config), linkage_config.dir];

return
end