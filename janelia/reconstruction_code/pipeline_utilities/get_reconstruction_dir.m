function reconstruction_dir = get_reconstruction_dir(config)
% image_dir = get_image_dir(config)
%
% Get the stack image directory from the config
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