function seg_param_choice = get_seg_param_choice_dir(config)
% image_dir = get_seg_param_choice_dir(config)
%
% Get the directory for segmentation parameter choices
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
%

stack_config = config.stack;
reconstruction_config = config.reconstruction;

curr_dir = pwd2;
cd(reconstruction_config.root_dir);
if(exist(stack_config.name, 'dir')~=7)
  mkdir2(stack_config.name)
  system(['chmod -R a+rw ', stack_config.name]); 
end;
cd(curr_dir);

reconstruction_dir = [reconstruction_config.root_dir, stack_config.name, '/'];

return
end