function config = process_config_structure(config)
% config = process_config_structure(config)
%
% Generate dependent configuration attributes based on user-
% specified configuration attributes.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  01232009  init code
%

global config_global code_dir

if(~isfield(config, 'is_verbose'))
  config.is_verbose = true;
end

% Create a temporary working directory. This is deleted at the end of
% pipeline's execution.
if(isunix)
  tmp_dir_root = '/tmp/';
else
  tmp_dir_root = '../tmp/';
end
is_created_tmp_dir = false;
while(~is_created_tmp_dir)
  tmp_dir_name = [tmp_dir_root, 'em_temp_wd_', num2str(round(rand(1)*10000000)), '/'];
  if(exist(tmp_dir_name, 'dir')==7)
    continue;
  end
  check_for_dir(tmp_dir_name);
  config_global.temp_dir = tmp_dir_name;
  config.temp_dir = tmp_dir_name;
  is_created_tmp_dir = true;
end

config.region.region_structure = get_region_structure(config);

config = initialize_region_config_information(config);

config.job.is_stand_alone = config_global.job.is_stand_alone;
config.code_dir = code_dir;
config.bin_dir = [code_dir, 'bin/'];

config_global.code_dir = config.code_dir;
config_global.bin_dir = config.bin_dir;

config_global.notify_email.gmailUser = 'em.reconstruct';
config_global.notify_email.gmailPassword = 'janelia.em.reconstruct';
config_global.notify_email.recipient = 'vitaladevunis@janelia.hhmi.org';

if(isfield(config, 'image_stack') && isfield(config.image_stack, 'filter_version') && ...
    ~iscell(config.image_stack.filter_version))
  config.image_stack.filter_version = {config.image_stack.filter_version};
end
return;
end
