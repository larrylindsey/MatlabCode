function pipeline_sab(config_file_name, module_ids)
% pipeline_sab(config_file_name, module_ids)
% Calls different modules of the pipeline based on module_id
%
% module_id
%   0-999.99        For reconstruction of Serial Section Electron
%                     Micrographs. Type help pipeline_serial_section.
%   999.99-1999.99  For reconstruction of Block Face and Focus Ion Beam
%                     Micrographs. Type help pipeline_block_face.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% AG mod for linkage
% v0  04162008  init code
% v1  01192009  pipeline version 3
% v2  03112009  branch into 2D and 3D reconstruction
% v3  01082010  stand alone binary
%

global config_global code_dir

l = load(config_file_name);
config = l.config;
if(config.is_verbose)
  fprintf('START: pipeline_sab\n');
  [s,r] = system('logname'); %#ok<ASGLU>
  fprintf('logname: %s\n', r(1:end-1));
  [s,r] = system('hostname'); %#ok<ASGLU>
  fprintf('hostname: %s\n', r(1:end-1));
end
config_global = l.config_global;
code_dir = l.code_dir;

config_global.temp_dir = './';
config.temp_dir = './';

if(~isfield(config, 'region') || ~isfield(config.region, 'dir'))
  config = initialize_region_config_information(config);
end

config.job.is_stand_alone = config_global.job.is_stand_alone;
config.code_dir = code_dir;
config.bin_dir = [code_dir, 'bin/'];

config_global.code_dir = config.code_dir;
config_global.bin_dir = config.bin_dir;

if(isfield(config, 'image_stack') && isfield(config.image_stack, 'filter_version') && ...
    ~iscell(config.image_stack.filter_version))
  config.image_stack.filter_version = {config.image_stack.filter_version};
end

fprintf('module_ids %s\n', module_ids);
module_ids = eval(module_ids);
module_ids = int32(module_ids*100);

for module_id = module_ids
  main_module_id = int32(floor(double(module_id)/100));
  sub_module_id = int32(mod(module_id,100));

  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Pipeline for Serial Section Electron Micrographs
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  if(main_module_id>=0 && main_module_id<999)
    pipeline_serial_section(config, main_module_id, sub_module_id);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Pipeline for FIB or Block Face Electron Micrographs
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  if(main_module_id>=1000 && main_module_id<2000)
    pipeline_block_face(config, main_module_id, sub_module_id);
  end
end;

if(config.is_verbose)
  fprintf('STOP: pipeline_sab\n');
end

return
end
