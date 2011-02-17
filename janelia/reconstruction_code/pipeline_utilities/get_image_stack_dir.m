function image_stack_dir = get_image_stack_dir(config)
% image_stack_dir = get_image_stack_dir(config)
% Get the directory for storing image sub-stacks
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  03112009  init code
%

image_stack_dir = [get_reconstruction_dir(config), config.image_stack.dir];
check_for_dir(image_stack_dir);
return
end
