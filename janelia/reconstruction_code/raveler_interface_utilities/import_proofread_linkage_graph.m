function import_proofread_linkage_graph(config)
% import_proofread_linkage_graph(config)
% Import linkage graph from proofread data
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  10222008  init code
%

stack_config = config.stack;
linkage_config = config.linkage;
import_config = config.proofreader.import;
raveler_config = import_config.Raveler;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('START: import_proofread_linkage_graph\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fin_sp_2_seg = fopen([raveler_config.proofread_data_dir, ...
  raveler_config.superpixel_to_segment_file_name], 'rt');
sp_2_seg = fscanf(fin_sp_2_seg, '%d', [3 inf])';
fclose(fin_sp_2_seg);

fin_seg_2_body = fopen([raveler_config.proofread_data_dir, ...
  raveler_config.segment_to_body_map_file_name], 'rt');
seg_2_body = fscanf(fin_seg_2_body, '%d', [2 inf])';
fclose(fin_seg_2_body);

%%%
% Find out the spatial extent along z-axis for all the bodies.
% If some body "skips" a section because it lies completely within a fold
% then create a dummy segment in the linkage graph to maintain continuity.
%%%
lists_of_bodies_create_dummy_segments_for = ...
  get_lists_of_bodies_to_create_dummy_segments_for(...
  uint32(sp_2_seg), uint32(seg_2_body));

seg_dir = [get_reconstruction_dir(config), ...
  config.superpixel(1).dir, import_config.dir];
seg_suffix = ['.prfrd_seg', '.', import_config.proofread_name];

image_prefixes_1 = {};
for i = 1:length(stack_config.case_ids)-1
  case_id = stack_config.case_ids(i);
  case_id_1 = stack_config.case_ids(i+1);
  fprintf('plane: %d\n', case_id);
  
  if(isempty(image_prefixes_1))
    image_prefixes_1 = get_image_prefixes_subdirs(config, case_id);
  end
  image_prefixes_2 = get_image_prefixes_subdirs(config, case_id_1);
  
  for tile_1 = 1:length(image_prefixes_1)
    file_name = [seg_dir, image_prefixes_1{tile_1}, ...
      seg_suffix, '.mat'];
    seg_1 = load2(file_name, 'superpixel_to_seg_label');
    seg_ids_1 = unique(nonzeros(seg_1.superpixel_to_seg_label(:,2)));
    seg_2_body_1 = seg_2_body(ismember(seg_2_body(:,1), seg_ids_1), :);
    seg_2_body_1 = sortrows(seg_2_body_1, 2);
    for tile_2 = 1:length(image_prefixes_2)
      fprintf('tile 1:%d, tile 2:%d\n', tile_1, tile_2);
      fprintf('%s\n%s\n', image_prefixes_1{tile_1}, image_prefixes_2{tile_2});
      file_name = [seg_dir, image_prefixes_2{tile_2}, ...
        seg_suffix, '.mat'];
      seg_2 = load2(file_name, 'superpixel_to_seg_label');
      seg_ids_2 = unique(nonzeros(seg_2.superpixel_to_seg_label(:,2)));
      seg_2_body_2 = seg_2_body(ismember(seg_2_body(:,1), seg_ids_2), :);
      seg_2_body_2 = sortrows(seg_2_body_2, 2);
      
      file_name_suffix = get_file_name_from_tuple(...
        [get_reconstruction_dir(config), linkage_config.dir], ...
        image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'lkp.');
      file_name = [file_name_suffix, '.proofread', ...
        '.', import_config.proofread_name, '.v0', ...
        config.segmentation_choose.choice.seg_suffix, '.mat'];
      if(exist(get_storage_file_name(file_name), 'file')==2)
        delete(get_storage_file_name(file_name));
      end
      
      if(isempty(seg_2_body_1) || isempty(seg_2_body_2))
        warning('Superpixel-to-segment map is empty!!');
        continue;
      end
      l_g = get_linkage_graph_from_seg_2_body_map(uint32(seg_2_body_1), ...
        uint32(seg_2_body_2));
      
      % add dummy segments for discontinuous bodies
      if(isfield(raveler_config, 'is_numbered_by_depth') && ...
          ~raveler_config.is_numbered_by_depth)
        is_disappearing = ismember(seg_2_body_1(:,2), ...
          lists_of_bodies_create_dummy_segments_for{config.region.case_ids==case_id_1});
      else
        is_disappearing = ismember(seg_2_body_1(:,2), ...
          lists_of_bodies_create_dummy_segments_for{case_id_1});
      end
      l_g_dummy_up = [seg_2_body_1(is_disappearing, 1), ...
        -seg_2_body_1(is_disappearing, 2)];
      if(isfield(raveler_config, 'is_numbered_by_depth') && ...
          ~raveler_config.is_numbered_by_depth)
        is_disappearing = ismember(seg_2_body_2(:,2), ...
          lists_of_bodies_create_dummy_segments_for{config.region.case_ids==case_id});
      else
        is_disappearing = ismember(seg_2_body_2(:,2), ...
          lists_of_bodies_create_dummy_segments_for{case_id});
      end
      l_g_dummy_down = [-seg_2_body_2(is_disappearing, 2), ...
        seg_2_body_2(is_disappearing, 1)];
      
      links_3D_p = [double(l_g), ones(size(l_g,1),1)*2^14; ...
        double(l_g_dummy_up), ones(size(l_g_dummy_up,1),1)*2^14;
        double(l_g_dummy_down), ones(size(l_g_dummy_down,1),1)*2^14];
      
      if(~isempty(links_3D_p))
        fprintf('Saving linkage graph:\n%s\n%s\n', file_name, ...
          get_storage_file_name(file_name));
        save2(file_name, 'links_3D_p');
      end
    end
  end
  
  image_prefixes_1 = image_prefixes_2;
end
fprintf('STOP: import_proofread_linkage_graph\n');
return
end
