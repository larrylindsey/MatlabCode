function set_segment_suffix_choice_to_proofread(config)
% set_segment_suffix_choice_to_proofread(config)
% Set superpixel choice to proofread suffix. This is useful for rest of the
% modules.
%
% Shiv N. Vitaladevuni
% v0  04062009  init. code
%

fprintf('START:set_segment_suffix_choice_to_proofread\n');

stack_config = config.stack;
import_config = config.proofreader.import;
choose_config = config.segmentation_choose;

superpixel_param_choice_dir = get_segmentation_param_choice_dir(config);

choice_param.method = 'proofread_result';
choice_param.seg_suffix = ['.prfrd_seg', '.', import_config.proofread_name];

for layer_id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(layer_id);
  fprintf('plane: %d\n', case_id);
  [image_prefixes, image_sub_dirs] = ...
    get_image_prefixes_subdirs(config, case_id);
  
  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))
      continue;
    end
    fprintf('tile: %d\n%s\n', tile_id, image_prefixes{tile_id});
    file_name = [superpixel_param_choice_dir, image_prefixes{tile_id}, ...
      choose_config.save_suffix, '.mat'];
    check_for_dir([superpixel_param_choice_dir, image_sub_dirs{tile_id}]);
    save2(file_name, 'choice_param');
  end
  fprintf('done.\n');
end

fprintf('STOP:set_segment_suffix_choice_to_proofread\n');
return
end
