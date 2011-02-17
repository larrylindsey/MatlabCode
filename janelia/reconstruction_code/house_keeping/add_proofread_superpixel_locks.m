function add_proofread_superpixel_locks(config)
% add_proofread_superpixel_locks(config)

fprintf('START: add_proofread_superpixel_locks\n');
stack_config = config.stack;
import_config = config.proofreader.import;

save_dir = [get_reconstruction_dir(config), ...
  config.superpixel(1).dir, import_config.dir];
save_suffix_c = ['.prfrd_csp', '.', import_config.proofread_name];
save_suffix_a = ['.prfrd_asp', '.', import_config.proofread_name];
for layer_id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(layer_id);
  fprintf('plane: %d\n', case_id);
  image_prefixes = get_image_prefixes_subdirs(config, case_id);
  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))
      continue;
    end
    fprintf('tile: %d\n%s\n', tile_id, image_prefixes{tile_id});

    file_name = [save_dir, image_prefixes{tile_id}, ...
      save_suffix_c, '.mat'];
    seg = load2(file_name);
    if(isfield(seg, 'locked_labels'))
      continue;
    end
    locked_labels = unique(seg.label_map(seg.label_map>0)); %#ok<NASGU>
    label_map = seg.label_map; %#ok<NASGU>
    save2(file_name, 'label_map', 'locked_labels');

    file_name = [save_dir, image_prefixes{tile_id}, ...
      save_suffix_a, '.mat'];
    seg = load2(file_name);
    if(isfield(seg, 'locked_labels'))
      continue;
    end
    locked_labels = unique(seg.label_map(seg.label_map>0)); %#ok<NASGU>
    label_map = seg.label_map; %#ok<NASGU>
    save2(file_name, 'label_map', 'locked_labels');
  end  
end

fprintf('STOP: add_proofread_superpixel_locks\n');
return
end
