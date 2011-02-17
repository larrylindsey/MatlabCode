function import_proofread_superpixel_map_uni_tile(config)
% import_proofread_superpixel_map_uni_tile(config)
% Import suerpixel maps for uni tile stacks. No need for alignment or
% adjusting.
% Label '0' is written on all undefined pixels.
%
% Shiv N. Vitaladevuni
% v0  09292009  init. code
%

fprintf('START: import_proofread_superpixel_map_uni_tile\n');
if(strcmp(config.proofreader.method, 'Raveler')~=1)
  error('Proofreader method is not Raveler');
end

fprintf('Importing proofread data from proofreader ..\n');

stack_config = config.stack;
import_config = config.proofreader.import;
raveler_config = import_config.Raveler;
if(~isfield(import_config, 'is_verbose'))
  import_config.is_verbose = true;
end
if(~isfield(import_config, 'is_verbose_figures'))
  import_config.is_verbose_figures = false;
end

save_dir = [get_reconstruction_dir(config), ...
  config.superpixel(1).dir, import_config.dir];
check_for_dir(save_dir);
save_suffix = ['.prfrd_asp', '.', import_config.proofread_name];
for layer_id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(layer_id);
  fprintf('plane: %d\n', case_id);
  [image_prefixes, image_sub_dirs] = ...
    get_image_prefixes_subdirs(config, case_id);
  image_prefix = image_prefixes{1};
  image_sub_dir = image_sub_dirs{1};
  fprintf('image_prefix: %s\n', image_prefix);

  label_map = imread(sprintf([raveler_config.proofread_data_dir, ...
    raveler_config.superpixel_dir, raveler_config.superpixel_prefix, ...
    raveler_config.superpixel_suffix], case_id));

  label_map(label_map<0) = 0;
  
  fprintf('Saving.\n');
  check_for_dir([save_dir, image_sub_dir]);
  file_name = [save_dir, image_prefix, save_suffix, '.mat'];
  locked_labels = unique(label_map(label_map>0)); %#ok<NASGU>
  save2(file_name, 'label_map', 'locked_labels');
  fprintf('done.\n');
end

fprintf('STOP:import_proofread_superpixel_map_uni_tile\n');
return
end
