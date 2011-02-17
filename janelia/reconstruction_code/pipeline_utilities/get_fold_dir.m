function fold_dir = get_fold_dir(config)
% fold_dir = get_fold_dir(config)
%
% Get the fold directory from the config
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  01022009  init code
%

fold_config = config.fold;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(fold_config.dir, 'dir')~=7)
  mkdir2(fold_config.dir);
end;
cd(prev_dir);

fold_dir = [get_reconstruction_dir(config), fold_config.dir];

return
end
