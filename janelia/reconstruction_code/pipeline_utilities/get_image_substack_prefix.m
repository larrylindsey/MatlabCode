function prefix = get_image_substack_prefix(config)
% prefix = get_image_substack_prefix(config)
%
% Get the prefix of file name of image sub-stack from the config
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
%

stack_config = config.stack;
image_stack_config = config.image_stack;

prefix = ['.', num2str(stack_config.case_ids(1)), '.', ...
  num2str(stack_config.case_ids(end)), ...
  sprintf('.%s', image_stack_config.filter_version{1:end})];
return
end