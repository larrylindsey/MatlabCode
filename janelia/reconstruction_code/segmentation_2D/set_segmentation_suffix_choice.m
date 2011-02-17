function set_segmentation_suffix_choice(config)
% set_segmentation_suffix_choice(config)(config)
% Set segmentation choice as specified in config.segmentation_choose.set_choice
% This is useful if blocks of sections have different parameters.
%
% Shiv N. Vitaladevuni
% v0  05212009  init. code
%

fprintf('START:set_segmentation_suffix_choice\n');

stack_config = config.stack;
choose_config = config.segmentation_choose;
choice_config = choose_config.set_choice;

segmentation_param_choice_dir = get_segmentation_param_choice_dir(config);

choice_param.method = choice_config.method;
choice_param.seg_suffix = choice_config.seg_suffix;

for layer_id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(layer_id);
  fprintf('plane: %d\n', case_id);
  [image_prefixes, image_sub_dirs, is_to_be_processed] = ...
    get_image_prefixes_subdirs(config, case_id);
  
  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))
      continue;
    end
    fprintf('tile: %d\n%s\n', tile_id, image_prefixes{tile_id});
    if(is_to_be_processed(tile_id)==0)
      fprintf('is_to_be_processed=false, skipping\n');
      continue;
    end
    
    file_name = [segmentation_param_choice_dir, image_prefixes{tile_id}, ...
      choose_config.save_suffix, '.mat'];
    check_for_dir([segmentation_param_choice_dir, image_sub_dirs{tile_id}]);
    save2(file_name, 'choice_param');
  end
  fprintf('done.\n');
end

fprintf('STOP:set_segmentation_suffix_choice\n');
return
end
