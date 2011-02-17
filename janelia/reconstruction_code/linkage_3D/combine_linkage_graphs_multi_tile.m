function combine_linkage_graphs_multi_tile(config)
% combine_linkage_graphs_multi_tile(config)
% For multi tile stacks, combine the linkage graphs of consecutive
% section pairs into one linkage graph.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  10222008  init code
%

if(config.is_verbose)
  fprintf('START: combine_linkage_graphs_multi_tile\n');
end
stack_config = config.stack;
linkage_config = config.linkage;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
apply_config = linkage_config.apply;

align_seg_dir = get_align_seg_dir(config);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Merging the linkage graphs from multiple tiles ..\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
label_pairs_overlap_t = cell(1, length(stack_config.case_ids)-1);
links_3D_t = cell(1, length(stack_config.case_ids)-1);
image_prefixes_1 = {};
for layer_id = 1:length(stack_config.case_ids)-1
  case_id = stack_config.case_ids(layer_id);
  case_id_1 = stack_config.case_ids(layer_id+1);
  fprintf('case_id: %d, case_id_1: %d\n', case_id, case_id_1);
  
  if(isempty(image_prefixes_1))
    [image_prefixes_1, image_sub_dirs_1] = get_image_prefixes_subdirs(config, case_id);
  end
  [image_prefixes_2, image_sub_dirs_2] = get_image_prefixes_subdirs(config, case_id_1);

  % load label mapping for case_id
  image_set_string = get_set_string(image_prefixes_1);
  save_dir = [align_seg_dir, image_sub_dirs_1{1}];
  label_mapping_0 = load2([save_dir, 'sec.', num2str(case_id), ...
    '.aligned_seg_mapping', image_set_string, ...
    config.segmentation_choose.choice.seg_suffix, '.mat'], 'label_mappings');
  for tile = 1:length(image_prefixes_1)
    label_mapping_0.label_mappings{tile}(label_mapping_0.label_mappings{tile}<0)=0;
  end
  
  % load label mapping for case_id_1
  image_set_string = get_set_string(image_prefixes_2);
  save_dir = [align_seg_dir, image_sub_dirs_2{1}];
  label_mapping_1 = load2([save_dir, 'sec.', num2str(case_id_1), ...
    '.aligned_seg_mapping', image_set_string, ...
    config.segmentation_choose.choice.seg_suffix, '.mat'], 'label_mappings');
  for tile = 1:length(image_prefixes_2)
    label_mapping_1.label_mappings{tile}(label_mapping_1.label_mappings{tile}<0)=0;
  end
  
  links_3D_t{layer_id} = [];
  label_pairs_overlap_t{layer_id} = [];
  for tile_1 = 1:length(image_prefixes_1)
   for tile_2 = 1:length(image_prefixes_2)
     fprintf('tile 1:%d, tile 2:%d\n', tile_1, tile_2);
     fprintf('tile 1: %s\ntile 2: %s\n', image_prefixes_1{tile_1}, ...
       image_prefixes_2{tile_2});
     if(isempty(label_mapping_0.label_mappings{tile_1}) || ...
         isempty(label_mapping_1.label_mappings{tile_2}))
       fprintf('Label mappings empty.\n');
       continue;
     end
     file_name_prefix = get_file_name_from_tuple(...
       [get_reconstruction_dir(config), linkage_config.dir], ...
       image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'lkp.');
     file_name = [file_name_prefix, ...
       '.', model_config.type, model_config.suffix, ...
       '.', feature_config.type, feature_config.suffix, ...
       apply_config.model_suffix, ...
       config.segmentation_choose.choice.seg_suffix, '.mat'];
     try
       fprintf('Attempting to load %s ...\n', file_name);
       linkage_graph = load2(file_name, 'links_3D_p', 'label_pairs_overlap');
     catch %#ok<CTCH>
       fprintf('failed.\nProceeding nevertheless.\n');
       continue;
     end
     fprintf('loaded.\n');
     if(isempty(linkage_graph.links_3D_p))
       fprintf('Linkage graph empty - skipping.\n');
       continue;
     end
     linkage_graph.links_3D_p = linkage_graph.links_3D_p(...
       linkage_graph.links_3D_p(:,1)<length(label_mapping_0.label_mappings{tile_1}) & ...
       linkage_graph.links_3D_p(:,2)<length(label_mapping_1.label_mappings{tile_2}), ...
       :);
     to_change = linkage_graph.links_3D_p(:,1)>0;
     linkage_graph.links_3D_p(to_change,1) = ...
       label_mapping_0.label_mappings{tile_1}(...
       linkage_graph.links_3D_p(to_change,1)+1);
     to_change = linkage_graph.links_3D_p(:,2)>0;
     linkage_graph.links_3D_p(to_change,2) = ...
       label_mapping_1.label_mappings{tile_2}(...
       linkage_graph.links_3D_p(to_change,2)+1);
     
     fprintf('Adding links_3D_p to links_3D_t.\n');
     links_3D_t{layer_id} = [links_3D_t{layer_id}; linkage_graph.links_3D_p];
     label_pairs_overlap_t{layer_id} = [label_pairs_overlap_t{layer_id}; ...
       linkage_graph.label_pairs_overlap];
   end
  end
  links_3D_t{layer_id} = ...
    unique(sortrows(links_3D_t{layer_id}), 'rows', 'last');
  label_pairs_overlap_t{layer_id} = ...
    unique(sortrows(label_pairs_overlap_t{layer_id}), 'rows', 'last');
  
  % Remove segments that have no superpixel members (possibly removed
  % during roi masking)
  if(~isempty(links_3D_t{layer_id}))
    keep_flag = links_3D_t{layer_id}(:,1)~=0 & links_3D_t{layer_id}(:,2)~=0;
    links_3D_t{layer_id} = links_3D_t{layer_id}(keep_flag,:);
  end
  if(~isempty(label_pairs_overlap_t{layer_id}))
    keep_flag = label_pairs_overlap_t{layer_id}(:,1)~=0 & ...
      label_pairs_overlap_t{layer_id}(:,2)~=0;
    label_pairs_overlap_t{layer_id} = label_pairs_overlap_t{layer_id}(keep_flag,:);
  end
  
  image_prefixes_1 = image_prefixes_2;
  image_sub_dirs_1 = image_sub_dirs_2;
end
save2([get_region_dir(config), 'links_3D_t', '.', ...
  model_config.type, '.', feature_config.type, apply_config.model_suffix, ...
  config.segmentation_choose.choice.seg_suffix, '.mat'], ...
  'links_3D_t', 'label_pairs_overlap_t');

if(config.is_verbose)
  fprintf('STOP: combine_linkage_graphs_multi_tile\n');
end
return
end
