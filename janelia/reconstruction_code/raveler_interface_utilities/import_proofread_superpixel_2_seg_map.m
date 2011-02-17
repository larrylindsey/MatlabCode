function import_proofread_superpixel_2_seg_map(config)
% import_proofread_superpixel_2_seg_map(config)
% Import proofread superpixel_2_seg_map from proofreader.
%
% Shiv N. Vitaladevuni
% v0  04062009  init. code
%

fprintf('START: import_proofread_superpixel_2_seg_map\n');
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
import_suffix = ['.prfrd_asp', '.', import_config.proofread_name];
save_suffix = ['.prfrd_seg', '.', import_config.proofread_name];

fin_sp_2_seg = fopen([raveler_config.proofread_data_dir, ...
  raveler_config.superpixel_to_segment_file_name], 'rt');
sp_2_seg = fscanf(fin_sp_2_seg, '%d', [3 inf])';
fclose(fin_sp_2_seg);

% Sometimes two segments in a section may be combined by the proofreader
% using 3D link. If a body is occuring in only one section, such
% combinations should be put in the superpixel_2_seg_label rather than
% linkage graph.
fin_seg_2_body = fopen([raveler_config.proofread_data_dir, ...
  raveler_config.segment_to_body_map_file_name], 'rt');
seg_2_body = fscanf(fin_seg_2_body, '%d', [2 inf])';
fclose(fin_seg_2_body);
sp_2_seg = double(get_sp_2_seg_remap_single_section_body_segment(...
  uint32(sp_2_seg), uint32(seg_2_body)));

for case_id = stack_config.case_ids
  layer_id = find(config.region.case_ids==case_id);
  fprintf('case_id: %d, layer_id: %d\n', case_id, layer_id);
  [image_prefixes, image_sub_dirs] = ...
    get_image_prefixes_subdirs(config, case_id);

  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))
      continue;
    end
    if(import_config.is_verbose)
      fprintf('tile:%d\n', tile_id);
    end
    file_name = [save_dir, image_prefixes{tile_id}, ...
      import_suffix, '.mat'];
    seg = load2(file_name, 'label_map');
    sp_ids = unique(seg.label_map(:));
    if(isfield(raveler_config, 'is_numbered_by_depth') && ...
        ~raveler_config.is_numbered_by_depth)
      sp_2_seg_l = sp_2_seg(sp_2_seg(:,1)==layer_id, 2:3);
    else
      sp_2_seg_l = sp_2_seg(sp_2_seg(:,1)==case_id, 2:3);
    end
    superpixel_to_seg_label = sp_2_seg_l(...
      ismember(sp_2_seg_l(:,1), sp_ids), :);
    file_name = [save_dir, image_prefixes{tile_id}, ...
      save_suffix, '.mat'];
    check_for_dir([save_dir, image_sub_dirs{tile_id}]);
    locked_labels = nonzeros(unique(superpixel_to_seg_label(:,2))); %#ok<NASGU>
    save2(file_name, 'superpixel_to_seg_label', 'locked_labels');
  end
  fprintf('done.\n');
end

fprintf('STOP: import_proofread_superpixel_2_seg_map\n');
return
end
