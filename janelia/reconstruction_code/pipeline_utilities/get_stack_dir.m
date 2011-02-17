function stack_dir = get_stack_dir(config)
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

stack_dir = [stack_config.dir, stack_config.name, '/'];

return
end