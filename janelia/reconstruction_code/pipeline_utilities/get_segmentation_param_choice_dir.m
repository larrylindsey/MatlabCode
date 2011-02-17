function segmentation_param_choice_dir = get_segmentation_param_choice_dir(config)
% segmentation_param_choice_dir = get_segmentation_param_choice_dir(config)
%
% Get the directory for segmentation parameter choices
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  12162008  init code
%

choose_config = config.segmentation_choose;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(choose_config.dir, 'dir')~=7)
  mkdir2(choose_config.dir);
end;
cd(prev_dir);

segmentation_param_choice_dir = [get_reconstruction_dir(config), choose_config.dir];

return
end