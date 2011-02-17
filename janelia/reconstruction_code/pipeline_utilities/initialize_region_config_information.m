function config = initialize_region_config_information(config, ...
  to_create_new_region_in_registry)
% config = initialize_region_config_information(config)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  01232009  init code
%

if(config.is_verbose)
  fprintf('START: initialize_region_config_information\n');
end
if(nargin<2)
  to_create_new_region_in_registry = true;
end

stack_config = config.stack;

if(~isfield(config, 'region') || ~isfield(config.region, 'case_ids'))
  config.region.case_ids = stack_config.case_ids;
end

if(isfield(stack_config, 'image_structure') && ~isempty(stack_config.image_structure))
  image_structure_name = get_file_name_from_full_path(stack_config.image_structure);
  image_structure_name = image_structure_name(1:end-4);
  image_structure_name_cleaned = strrep(stack_config.image_structure, '/', '_');
  config.region.image_structure = stack_config.image_structure;
else
  image_structure_name = 'fs';
  image_structure_name_cleaned = 'fs';
  config.region.image_structure = '';
end
case_ids_string = [sprintf('%g_', config.region.case_ids(1:end-1)), ...
  num2str(config.region.case_ids(end), '%g')];

config.region.name_original =  ['region.', image_structure_name, ...
  '.', num2str(config.region.case_ids(1)), '.', num2str(config.region.case_ids(end)), ...
  '.', image_structure_name_cleaned, '.', case_ids_string];


region_registry_dir = get_region_registry_dir(config);
save_file_name_orig = [region_registry_dir, config.region.name_original, '.mat'];
if(to_create_new_region_in_registry)
  region_info = config.region; %#ok<NASGU>
  save_file_name_hashed = save2(save_file_name_orig, 'region_info');
  save_file_name_hashed = get_file_name_from_full_path(save_file_name_hashed);
  
  config.region.dir = [save_file_name_hashed(1:end-4), '/'];
  region_info = config.region; %#ok<NASGU>
  save_file_name_hashed_2 = save2(save_file_name_orig, 'region_info');
  save_file_name_hashed_2 = get_file_name_from_full_path(save_file_name_hashed_2);
  if(strcmp(save_file_name_hashed, save_file_name_hashed_2)~=1)
    error('Error while creating region directory: hash values are inconsistent.');
  end
  
  if(exist([get_reconstruction_dir(config), config.region.dir], 'dir')~=7)
    mkdir2([get_reconstruction_dir(config), config.region.dir]);
  end
  
  if(isfield(stack_config, 'image_structure') && ~isempty(stack_config.image_structure))
    [status, result] = system(['cp ', get_stack_dir(config), stack_config.image_structure, ' ', ...
      get_reconstruction_dir(config), config.region.dir, ...
      get_file_name_from_full_path(stack_config.image_structure)]);
    if(status~=0)
      fprintf(result);
      warning('Could not copy stack structure file to reconstruction region dir.');
    end
  end
else
  load2(save_file_name_orig, 'region_info');
  config.region.dir = region_info.dir; %#ok<NODEF>
end

if(config.is_verbose)
  fprintf('STOP: initialize_region_config_information\n');
end
return;
end
