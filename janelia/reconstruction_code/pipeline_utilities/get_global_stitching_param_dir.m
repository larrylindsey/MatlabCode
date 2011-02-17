function global_stitching_param_dir = get_global_stitching_param_dir(config)
% global_stitching_param_dir = get_global_stitching_param_dir(config)
%
% Get directory for storing global stitching parameters for a region.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
%

global_stitching_param_dir = [get_region_dir(config), ...
  'global_stitching_parameters/'];

return
end
