function region_registry_dir = get_region_registry_dir(config)
% segmentation_dir = get_region_registry_dir(config)
%
% Get the region registry directory from the config
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  01232009  init code
%

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(config.region.registry_dir, 'dir')~=7)
  mkdir2(config.region.registry_dir);
end;
cd(prev_dir);

region_registry_dir = ...
  [get_reconstruction_dir(config), config.region.registry_dir];

return
end