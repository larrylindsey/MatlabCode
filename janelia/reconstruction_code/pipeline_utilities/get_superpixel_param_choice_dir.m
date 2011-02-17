function superpixel_param_choice_dir = get_superpixel_param_choice_dir(config)
% image_dir = get_superpixel_param_choice_dir(config)
%
% Get the directory for superpixel parameter choices
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
%

choose_config = config.superpixel_choose;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(choose_config.dir, 'dir')~=7)
  mkdir2(choose_config.dir);
end;
cd(prev_dir);

superpixel_param_choice_dir = [get_reconstruction_dir(config), choose_config.dir];

return
end