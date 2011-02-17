function segmentation_dir = get_train_segmentation_dir(config)
% segmentation_dir = get_train_segmentation_dir(config)
%
% Get the segemntation training directory from the config
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
%

train_config = config.train_segmentation;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(train_config.dir, 'dir')~=7)
  mkdir2(train_config.dir);
end;
cd(prev_dir);

segmentation_dir = [get_reconstruction_dir(config), train_config.dir];

return
end