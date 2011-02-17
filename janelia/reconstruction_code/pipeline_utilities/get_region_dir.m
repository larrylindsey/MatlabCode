function region_dir = get_region_dir(config)
% region_dir = get_region_dir(config)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  01232009  init code
%

region_dir = [get_reconstruction_dir(config), config.region.dir];

return;
end
