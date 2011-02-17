function process_superpixel_2_2Dseg_linkage_graph(config)
% process_superpixel_2_2Dseg_linkage_graph(config)
%
% Makes sure the segment ids are unique for the whole stack
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
% 
% v0  ~03202008 init code
% v1  04142008  modification to reconstruction pipeline
% 

if(config.is_verbose)
  fprintf('START: process_superpixel_2_2Dseg_linkage_graph\n');
end

linkage_config = config.linkage;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
apply_config = linkage_config.apply;
stack_config = config.stack;

fprintf('Loading sectionwise superpixel-to-segment maps\n');
fprintf('Reading in superpixel-to-segment maps of sections ..\n');
superpixel_2_seg_map = {};
superpixel_ids_sets = {};
for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  sp_2_seg = ...
    load2([get_region_dir(config), 'superpixel_2_seg_map_t.', num2str(case_id), ...
    config.segmentation_choose.choice.seg_suffix, '.mat'], ...
    'superpixel_2_seg_map_t', 'superpixel_ids_sets');
  superpixel_2_seg_map{i} = sp_2_seg.superpixel_2_seg_map_t; %#ok<AGROW>
  fprintf('section: %d, min segment id: %u, max segment id: %u\n', ...
    case_id, min(superpixel_2_seg_map{i}), max(superpixel_2_seg_map{i}));
  superpixel_ids_sets{i} = sp_2_seg.superpixel_ids_sets; %#ok<AGROW>
end

fprintf('Loading preliminary linkage graph\n');
load2([get_region_dir(config), 'links_3D_t', '.', ...
  model_config.type, '.', feature_config.type, apply_config.model_suffix, ...
  config.segmentation_choose.choice.seg_suffix, '.mat'], 'links_3D_t', ...
  'label_pairs_overlap_t');
links_3D = links_3D_t;
label_pairs_overlap = label_pairs_overlap_t;
clear links_3D_t label_pairs_overlap_t

fprintf(['Removing links between segments that are no longer referenced by\n', ...
  'superpixel_2_seg_map\n']);
if(isempty(links_3D))
  links_3D = cell(1, length(superpixel_2_seg_map)-1);
  label_pairs_overlap = cell(1, length(superpixel_2_seg_map)-1);
else
  fprintf('Section: 1\n');
  if(~isempty(links_3D{1}))
    to_retain = ismember(links_3D{1}(:,1), superpixel_2_seg_map{1}) | ...
      links_3D{1}(:,1)<0;
    links_3D{1} = links_3D{1}(to_retain, :);
  end
  if(~isempty(label_pairs_overlap{1}))
    to_retain = ismember(label_pairs_overlap{1}(:,1), superpixel_2_seg_map{1});
    label_pairs_overlap{1} = label_pairs_overlap{1}(to_retain, :);
  end
  for i = 2:length(superpixel_2_seg_map)-1
    fprintf('Section : %d\n', i);
    if(~isempty(links_3D{i-1}))
      to_retain = ismember(links_3D{i-1}(:,2), superpixel_2_seg_map{i}) | ...
        links_3D{i-1}(:,2)<0;
      links_3D{i-1} = links_3D{i-1}(to_retain, :);
    end
    if(~isempty(links_3D{i}))
      to_retain = ismember(links_3D{i}(:,1), superpixel_2_seg_map{i}) | ...
        links_3D{i}(:,1)<0;
      links_3D{i} = links_3D{i}(to_retain, :);
    end
    if(~isempty(label_pairs_overlap{i-1}))
      to_retain = ismember(label_pairs_overlap{i-1}(:,2), superpixel_2_seg_map{i});
      label_pairs_overlap{i-1} = label_pairs_overlap{i-1}(to_retain, :);
    end
    if(~isempty(label_pairs_overlap{i}))
      to_retain = ismember(label_pairs_overlap{i}(:,1), superpixel_2_seg_map{i});
      label_pairs_overlap{i} = label_pairs_overlap{i}(to_retain, :);
    end
  end
  fprintf('Section: %d\n', length(superpixel_2_seg_map));
  if(~isempty(links_3D{end}))
    to_retain = ismember(links_3D{end}(:,2), superpixel_2_seg_map{end}) | ...
      links_3D{end}(:,2)<0;
    links_3D{end} = links_3D{end}(to_retain, :);
  end
  if(~isempty(label_pairs_overlap{end}))
    to_retain = ismember(label_pairs_overlap{end}(:,2), superpixel_2_seg_map{end});
    label_pairs_overlap{end} = label_pairs_overlap{end}(to_retain, :);
  end
end

fprintf('Relabeling segments to remove gaps in ids ..\n');
annotations_voxel = {};
annotations_body_segment = {};
fprintf('case_id: %d\n', stack_config.case_ids(1));
seg_ids = unique(nonzeros(superpixel_2_seg_map{1}));
seg_relabel = [];
seg_relabel(1) = 0;
seg_relabel(seg_ids+1) = 1:length(seg_ids);
superpixel_2_seg_map{1} = seg_relabel(1+superpixel_2_seg_map{1});
if(~isempty(links_3D{1}))
  to_change = links_3D{1}(:,1)>0;
  links_3D{1}(to_change,1) = seg_relabel(1+links_3D{1}(to_change,1));
end
if(~isempty(label_pairs_overlap{1}))
  label_pairs_overlap{1}(:,1) = seg_relabel(1+label_pairs_overlap{1}(:,1));
end
annotations_file_name = [get_region_dir(config), config.annotations.dir, ...
  'annotations.', num2str(stack_config.case_ids(1)), '.txt'];
if(exist(annotations_file_name, 'file')==2)
  fprintf('reading: %s\n', annotations_file_name);
  [annotations_voxel{1}, annotations_body_segment{1}] = ...
    get_tile_voxel_annotations(annotations_file_name);
  if(~isempty(annotations_voxel{1}))
    for j = 1:length(annotations_voxel{1})
      annotations_voxel{1}(j).z = stack_config.case_ids(1); %#ok<AGROW>
      annotations_voxel{1}(j).annotation = transform_tbar_postsynaptic_coordinates(...
        annotations_voxel{1}(j).annotation, [1, 0, 0; 0, 1, 0], [], [], stack_config.case_ids(1)); %#ok<AGROW>
    end
  end
  if(~isempty(annotations_body_segment{1}))
    for j = 1:length(annotations_body_segment{1})
      annotations_body_segment{1}(j).segments = seg_relabel(...
        1 + annotations_body_segment{1}(j).segments); %#ok<AGROW>
    end
  end
end
for i = 2:length(superpixel_2_seg_map)-1
  case_id = stack_config.case_ids(i);
  fprintf('case_id: %d\n', case_id);
  
  seg_ids = unique(nonzeros(superpixel_2_seg_map{i}));
  seg_relabel = [];
  seg_relabel(1) = 0;
  seg_relabel(seg_ids+1) = 1:length(seg_ids); %#ok<AGROW>
  seg_relabel(1) = 0;
  superpixel_2_seg_map{i} = seg_relabel(1+superpixel_2_seg_map{i}); %#ok<AGROW>

  if(~isempty(links_3D{i-1}))
    to_change = links_3D{i-1}(:,2)>0;
    links_3D{i-1}(to_change,2) = seg_relabel(1+links_3D{i-1}(to_change,2));
  end
  if(~isempty(links_3D{i}))
    to_change = links_3D{i}(:,1)>0;
    links_3D{i}(to_change,1) = seg_relabel(1+links_3D{i}(to_change,1));
  end
  if(~isempty(label_pairs_overlap{i-1}))
    label_pairs_overlap{i-1}(:,2) = seg_relabel(1+label_pairs_overlap{i-1}(:,2));
  end
  if(~isempty(label_pairs_overlap{i}))
    label_pairs_overlap{i}(:,1) = seg_relabel(1+label_pairs_overlap{i}(:,1));
  end
  
  annotations_file_name = [get_region_dir(config), config.annotations.dir, ...
    'annotations.', num2str(case_id), '.txt'];
  if(exist(annotations_file_name, 'file')==2)
    fprintf('reading: %s\n', annotations_file_name);
    [annotations_voxel{i}, annotations_body_segment{i}] = ...
      get_tile_voxel_annotations(annotations_file_name); %#ok<AGROW>
    if(~isempty(annotations_voxel{i}))
      for j = 1:length(annotations_voxel{i})
        annotations_voxel{i}(j).z = case_id; %#ok<AGROW>
        annotations_voxel{i}(j).annotation = transform_tbar_postsynaptic_coordinates(...
          annotations_voxel{i}(j).annotation, [1, 0, 0; 0, 1, 0], [], [], case_id); %#ok<AGROW>
      end
    end
    if(~isempty(annotations_body_segment{i}))
      for j = 1:length(annotations_body_segment{i})
        annotations_body_segment{i}(j).segments = seg_relabel(...
          1 + annotations_body_segment{i}(j).segments); %#ok<AGROW>
      end
    end
  end
end
i = length(superpixel_2_seg_map);
case_id = stack_config.case_ids(i);
fprintf('case_id: %d\n', case_id);
seg_ids = unique(nonzeros(superpixel_2_seg_map{end}));
seg_relabel = [];
seg_relabel(1) = 0;
seg_relabel(seg_ids+1) = 1:length(seg_ids);
superpixel_2_seg_map{end} = seg_relabel(1+superpixel_2_seg_map{end});
if(~isempty(links_3D{end}))
  to_change = links_3D{end}(:,2)>0;
  links_3D{end}(to_change,2) = seg_relabel(1+links_3D{end}(to_change,2));
end
if(~isempty(label_pairs_overlap{end}))
  label_pairs_overlap{end}(:,2) = seg_relabel(1+label_pairs_overlap{end}(:,2));
end
annotations_file_name = [get_region_dir(config), config.annotations.dir, ...
  'annotations.', num2str(case_id), '.txt'];
if(exist(annotations_file_name, 'file')==2)
  fprintf('reading: %s\n', annotations_file_name);
  [annotations_voxel{i}, annotations_body_segment{i}] = ...
    get_tile_voxel_annotations(annotations_file_name);
  if(~isempty(annotations_voxel{i}))
    for j = 1:length(annotations_voxel{i})
      annotations_voxel{i}(j).z = case_id; %#ok<AGROW>
      annotations_voxel{i}(j).annotation = transform_tbar_postsynaptic_coordinates(...
        annotations_voxel{i}(j).annotation, [1, 0, 0; 0, 1, 0], [], [], case_id); %#ok<AGROW>
    end
  end
  if(~isempty(annotations_body_segment{i}))
    for j = 1:length(annotations_body_segment{i})
      annotations_body_segment{i}(j).segments = seg_relabel(...
        1 + annotations_body_segment{i}(j).segments); %#ok<AGROW>
    end
  end
end

fprintf(['Making sure that each segment id is unique for the entire stack\n', ...
  'this would affect links_3D and annotations too.\n']);
seg_id_offset = 0;
for i = 2:length(superpixel_2_seg_map)-1
  fprintf('Section: %d\n', i);
  case_id = stack_config.case_ids(i);
  if(~isempty(superpixel_2_seg_map{i-1}))
    seg_id_offset = max(seg_id_offset, max(superpixel_2_seg_map{i-1})) + 1;
  else
    seg_id_offset = seg_id_offset + 1;
  end
  
  if(~isempty(superpixel_2_seg_map{i}))
    superpixel_2_seg_map{i}(superpixel_2_seg_map{i}>0) = ...
      superpixel_2_seg_map{i}(superpixel_2_seg_map{i}>0) + seg_id_offset; %#ok<AGROW> % don't increment for 0 label
  end
  fprintf('section: %d, min segment id: %u, max segment id: %u\n', ...
    case_id, min(superpixel_2_seg_map{i}), max(superpixel_2_seg_map{i}));
  
  if(~isempty(links_3D{i-1}))
    to_change = links_3D{i-1}(:,2)>0;
    links_3D{i-1}(to_change,2) = links_3D{i-1}(to_change,2) + seg_id_offset;
  end
  if(~isempty(label_pairs_overlap{i-1}))
    label_pairs_overlap{i-1}(:,2) = label_pairs_overlap{i-1}(:,2) + seg_id_offset;
  end
  if(~isempty(links_3D{i}))
    to_change = links_3D{i}(:,1)>0;
    links_3D{i}(to_change,1) = links_3D{i}(to_change,1) + seg_id_offset;
  end
  if(~isempty(label_pairs_overlap{i}))
    label_pairs_overlap{i}(:,1) = label_pairs_overlap{i}(:,1) + seg_id_offset;
  end
  
  if(~isempty(annotations_body_segment) && ~isempty(annotations_body_segment{i}))
    for j = 1:length(annotations_body_segment{i})
      annotations_body_segment{i}(j).segments = seg_id_offset + ...
        annotations_body_segment{i}(j).segments; %#ok<AGROW>
    end
  end
end;
if(length(superpixel_2_seg_map)>1)
  seg_id_offset = max(seg_id_offset, max(superpixel_2_seg_map{end-1})) + 1;
  superpixel_2_seg_map{end}(superpixel_2_seg_map{end}>0) = ...
    superpixel_2_seg_map{end}(superpixel_2_seg_map{end}>0) + seg_id_offset; % don't increment for 0 label
  if(~isempty(links_3D{end}))
    to_change = links_3D{end}(:,2)>0;
    links_3D{end}(to_change,2) = links_3D{end}(to_change,2) + seg_id_offset;
  end
  if(~isempty(label_pairs_overlap{end}))
    label_pairs_overlap{end}(:,2) = label_pairs_overlap{end}(:,2) + seg_id_offset;
  end
  if(~isempty(annotations_body_segment) && ~isempty(annotations_body_segment{end}))
    for j = 1:length(annotations_body_segment{end})
      annotations_body_segment{end}(j).segments = seg_id_offset + ...
        annotations_body_segment{end}(j).segments; %#ok<AGROW>
    end
  end
end

fprintf('Saving superpixel-to-segment map for Raveler\n');
save2([get_region_dir(config), 'superpixel_2_seg_map', ...
  config.segmentation_choose.choice.seg_suffix, '.mat'], ...
  'superpixel_2_seg_map');
output_to_proofreader_superpixel_2_seg_map(superpixel_2_seg_map, ...
  superpixel_ids_sets, config);

fprintf('Saving segment-to-body map for Raveler\n');
save2([get_region_dir(config), 'links_3D',  '.',  model_config.type, '.', ...
  feature_config.type, apply_config.model_suffix, ...
  config.segmentation_choose.choice.seg_suffix, '.mat'], ...
  'links_3D', 'label_pairs_overlap');

seg_2_body_map = output_to_proofreader_linkage_graph(...
  links_3D, superpixel_2_seg_map, label_pairs_overlap, config);

% map body-segment annotations to body and take their union
annotations_body = [];
if(~isempty(annotations_body_segment))
  for i = 1:length(annotations_body_segment)
    for j = 1:length(annotations_body_segment{i})
      b_ids = unique(seg_2_body_map(1+annotations_body_segment{i}(j).segments));
      for k = b_ids
        if(isempty(annotations_body))
          annotations_body.body_id = k;
          annotations_body.annotation{1} = ...
            annotations_body_segment{i}(j).annotation;
        else
          a_id = find([annotations_body(:).body_id]==k, 1);
          if(~isempty(a_id))
            annotations_body(a_id).annotation = union(...
              annotations_body(a_id).annotation, ...
              {annotations_body_segment{i}(j).annotation}); %#ok<AGROW>
          else
            annotations_body(end+1).body_id = k; %#ok<AGROW>
            annotations_body(end).annotation{1} = ...
              annotations_body_segment{i}(j).annotation; %#ok<AGROW>
          end
        end
      end
    end
  end
end

annotations_file_name = [get_to_be_proofread_dir(config), 'annotations.txt'];
output_stack_annotations(annotations_voxel, ...
  annotations_body, annotations_file_name, config);

if(copyfile(...
    [get_reconstruction_dir(config), linkage_config.dir, 'link_model', '.', ...
    model_config.type, '.', feature_config.type, apply_config.model_suffix, '.mat'], ...
    [get_to_be_proofread_dir(config), 'link_model', '.', model_config.type, '.', ...
    feature_config.type, apply_config.model_suffix, '.mat'],'f')==0)
  warning('Could not copy linkage model file'); %#ok<WNTAG>
end
if(config.is_verbose)
  fprintf('STOP: process_superpixel_2_2Dseg_linkage_graph\n');
end
return
end
