function gen_link_gph_overlap_hist_boundary_intr_hist_lp_wekarf_dmesh(config)
% gen_link_gph_overlap_hist_boundary_intr_hist_lp_dmesh(config)
% generate linkage graph using boosted classifier trained on
% overlap hist, intensity histograms features
% Section alignment : Deformable mesh - see
% align_stack_deformable_mesh_tile_pair_inter_plane.m
%
% - Takes as input a stack of images and 2D superpixel and segment maps.
% - Uses trained classifier from linkage_3D_train_intensity_pair_boost.m to perform 3D linkage.
% - Constructs a linkage graph, links_3D{}, with segments as nodes and 3D link
% confidences as the edge weights. This graph's link-weights are
% thresholded in proofread_linkage_superpixel.m to get the final 3D
% segmentation.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

global config_global

stack_config = config.stack;
linkage_config = config.linkage;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
apply_config = linkage_config.apply;
dmesh_config = config.align.linkage_align.deformable_mesh;

if(linkage_config.is_verbose)
  fprintf('START: gen_link_gph_overlap_hist_boundary_intr_hist_lp_wekarf_dmesh\n');
end

version_names{1} = 'overlap_hist_boundary_interior_hist';
version_names{2} = 'overlap_hist_boundary_interior_hist_shrink_LS';
version_id = find(strcmp(feature_config.type, version_names)); %#ok<EFIND>
if(isempty(version_id))
  error('Feature type does not match with called function. Exiting');
end;
if(strcmp(model_config.type, 'wekarf_lp_p3')==0)
  error('Classifier type does not match with called function. Exiting');
end;

get_linkage_dir(config);

dmesh_dir = get_deformable_mesh_dir(config);

junk = tree_node_w(1); %#ok<NASGU>

linkage_model = load2([get_reconstruction_dir(config), linkage_config.dir, 'link_model', '.', ...
  'boost_lp_p3', model_config.suffix, '.', ...
  feature_config.type, feature_config.suffix, ...
  apply_config.model_suffix, '.mat'], ...
  'classifier', 'n_dimension', 'classifier_insection', 'n_dimension_insection', ...
  'model_config', 'feature_config');

linkage_model_wekarf = load('/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/scripts/linkage_experiments/classification_exp/linkage_weka_random_forest.596.598.mat');

weak_learner = tree_node_w(linkage_model.model_config.tree_depth); %#ok<NASGU>

case_id = stack_config.case_ids(1);
fprintf('Constructing linkage graph (links_3D) ..\n%d ', case_id);

[images_1, image_prefixes_1, image_sub_dirs_1, is_to_be_processed_1, size_images_1] = ...
  get_image_from_stack(config, case_id); %#ok<ASGLU>
if(isfield(feature_config, 'filter_version') && ...
    ~isempty(feature_config.filter_version))
  for tile = 5%1:length(images_1)
    images_1{tile} = filter_image2(images_1{tile}, feature_config.filter_version, config);
  end
end
segment_2D_label_map_1 = load_segment_maps(image_prefixes_1, size_images_1, config);

for i = 2:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  fprintf('case_id: %d\n', case_id);

  [images_2, image_prefixes_2, image_sub_dirs_2, ...
    is_to_be_processed_2, size_images_2] = get_image_from_stack(config, case_id); %#ok<ASGLU>
  
  if(isfield(feature_config, 'filter_version') && ...
      ~isempty(feature_config.filter_version))
    for tile = 9%1:length(image_prefixes_2)
      images_2{tile} = filter_image2(images_2{tile}, feature_config.filter_version, config);
    end
  end
  segment_2D_label_map_2 = load_segment_maps(image_prefixes_2, size_images_2, config);

  % compute a linkage graph for each pair of overlapping tiles
  for tile_1 = 5%1:length(image_prefixes_1)
    for tile_2 = 9%1:length(image_prefixes_2)
      fprintf('--- Tile pair ---\n%d,%d\n', tile_1, tile_2);
      fprintf('%s\n%s\n', image_prefixes_1{tile_1}, image_prefixes_2{tile_2});
      
      if(is_to_be_processed_1(tile_1)==0 && is_to_be_processed_2(tile_2)==0)
        fprintf('Both tiles have is_to_be_processed=false, skipping this pair\n');
        continue;
      end
      
      file_name_prefix = get_file_name_from_tuple(...
        [get_reconstruction_dir(config), linkage_config.dir], ...
        image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'lkp.');
      file_name_suffix = ['.', model_config.type, model_config.suffix, ...
        '.', feature_config.type, feature_config.suffix, ...
        apply_config.model_suffix];
      file_name_linkage_graph = [file_name_prefix, file_name_suffix, ...
        config.segmentation_choose.choice.seg_suffix, '.mat'];
      fprintf('Deleting the target linkage file if it already exists\n');
      delete(get_storage_file_name(file_name_linkage_graph));
      
      if(isfield(apply_config, 'delete_linkage_graph_dont_create') && ...
          apply_config.delete_linkage_graph_dont_create)
        continue;
      end
      
      fprintf('Loading transforms between tile 1 and tile 2\n');
      [transforms_tp, transforms_tp_rev] = ...
        load_tile_pair_deformable_mesh_transforms(image_prefixes_1{tile_1}, ...
        image_prefixes_2{tile_2}, dmesh_config, dmesh_dir);
      
      links_temp = [];
      
      fprintf('Relabeling to remove gaps in label indexes and make non-overlapping\n');
      [seg_1_original, relabeling_1] = relabel_to_remove_nonexistent_labels(...
        segment_2D_label_map_1{tile_1}.original_label_map);
      [seg_2_original, relabeling_2] = relabel_to_remove_nonexistent_labels(...
        segment_2D_label_map_2{tile_2}.original_label_map);
      seg_2_original(seg_2_original>0) = seg_2_original(seg_2_original>0) + ...
        max(seg_1_original(:));
      relabeling_2(relabeling_2>0) = relabeling_2(relabeling_2>0) + ...
        max(seg_1_original(:));
      
      seg_1 = relabeling_1(1+segment_2D_label_map_1{tile_1}.label_map);
      seg_2 = relabeling_2(1+segment_2D_label_map_2{tile_2}.label_map);
      
      label_backref_1 = [];
      label_backref_1(1+relabeling_1(relabeling_1>0)) = ...
        find(relabeling_1)-1; %#ok<AGROW>
      label_backref_2 = [];
      label_backref_2(1+relabeling_2(relabeling_2>0)) = ...
        find(relabeling_2)-1; %#ok<AGROW>
      if(nnz(label_backref_1(1+seg_1_original)~=...
          segment_2D_label_map_1{tile_1}.original_label_map)~=0)
        error('Error in relabeling and backref of tile 1');
      end
      if(nnz(label_backref_2(1+seg_2_original)~=...
          segment_2D_label_map_2{tile_2}.original_label_map)~=0)
        error('Error in relabeling and backref of tile 2');
      end
      
      fprintf('Computing linkage confidences\n');
      if(isfield(linkage_config, 'modified_segment') || ...
          isfield(apply_config, 'retain_links_subset'))
        [features, label_pairs, label_pairs_original] = ...
          get_linkage_features_overlap_hist_shrink_LS(...
          seg_1, images_1{tile_1}, seg_2, images_2{tile_2}, ...
          transforms_tp, transforms_tp_rev, feature_config, ...
          seg_1_original, seg_2_original);
        if(~isempty(label_pairs_original))
          label_pairs_original = label_pairs_original(...
            label_pairs_original(:,1)>0 & label_pairs_original(:,2)>0, :);
        end
      else
        [features, label_pairs] = get_linkage_features_overlap_hist_shrink_LS(...
          seg_1, images_1{tile_1}, seg_2, images_2{tile_2}, ...
          transforms_tp, transforms_tp_rev, feature_config);
      end
      if(~isempty(features))
        write_arff('link_temp', features, ones(size(features,1),1), ...
          [config_global.temp_dir, 'link_temp.arff']);
        data_weka = loadARFF([config_global.temp_dir, 'link_temp.arff']);
        [junk, link_confidence] = wekaClassify(data_weka, linkage_model_wekarf.classifier_random_forest); %#ok<ASGLU>
        links_temp = [label_pairs, link_confidence(:,1) - 0.5];
        links_temp = links_temp(links_temp(:,3)>0.3 | links_temp(:,3)<-0.05, :);
        links_temp(:,3) = links_temp(:,3)*10;
      end

      if(~isempty(links_temp))
        links_temp = sortrows(links_temp);
        
        if(isfield(apply_config, 'retain_links_subset'))
          switch(apply_config.retain_links_subset)
            case 'max_confidence'
              l1 = links_temp(links_temp(:,3)>0, :);
              if(~isempty(l1))
                l1 = sortrows(l1, 3);
                [junk, I] = unique(l1(:,1), 'last'); %#ok<ASGLU>
                l1 = l1(I,:);
              end
              
              l2 = links_temp(links_temp(:,3)>0, :);
              if(~isempty(l2))
                l2 = sortrows(l2, 3);
                [junk, I] = unique(l2(:,2), 'last'); %#ok<ASGLU>
                l2 = l2(I,:);
              end
              
              l3 = links_temp(links_temp(:,3)>2.5, :);
              
              l4 = links_temp(links_temp(:,3)<=0, :);
              
              links_temp = [l1; l2; l3; l4];
              links_temp = unique(links_temp, 'rows');
          end
        end
      end
      
      [features, label_pairs] = get_linkage_features_overlap_area_shrink_LS(...
        seg_1, images_1{tile_1}, seg_2, images_2{tile_2}, ...
        transforms_tp, transforms_tp_rev);
      if(~isempty(features))
        link_confidence = (features(:,1)>10000).*(1 + features(:,1)/10000);
        links_temp = [links_temp; label_pairs, link_confidence]; %#ok<AGROW>
      end
      links_boost_area = links_temp;
      
      if(~isempty(links_temp))
        links_temp = links_temp(links_temp(:,1)>0 & links_temp(:,2)>0, :);
        
        % For each pair of segments, retain the confidence with maximum
        % absolute value.
        link_conf_abs = abs(links_temp(:,3));
        [junk, I] = sort(link_conf_abs); %#ok<ASGLU>
        links_temp = links_temp(I,:);
        [junk, I] = unique(links_temp(:, [1 2]), 'rows', 'last'); %#ok<ASGLU>
        links_temp = links_temp(I,:);
      end
      
      if(~isempty(links_temp))
        link_positive_threshold = 2;
        links_temp = links_temp(links_temp(:,3)>link_positive_threshold | ...
          links_temp(:,3)<-1, :);
        
        fprintf('Compute within section merger confidences\n');
        boundary_hist_1 = collect_segment_pair_stats_boundary_hist(...
          uint32(seg_1_original), uint8(255*images_1{tile_1}), ...
          uint8(linkage_model.feature_config.boundary_hist.bin_size));
        boundary_hist_1 = boundary_hist_1(min(boundary_hist_1(:,1:2), [], 2)>0, :);
        boundary_hist_1 = boundary_hist_1(...
          boundary_hist_1(:,1)~=boundary_hist_1(:,2), :);
        features_boundary = interval_sum(boundary_hist_1(:, 3:end));
        
        interior_hist_t = collect_segment_stats_interior_hist(...
          uint32(seg_1_original), uint8(255*images_1{tile_1}), ...
          uint8(feature_config.interior_hist.bin_size));
        interior_hist = [];
        interior_hist(interior_hist_t(:,1)+1, :) = interior_hist_t(:, 2:end); %#ok<AGROW>
        interior_hist = interval_sum(interior_hist);
        features_interior = [interior_hist(boundary_hist_1(:,1)+1, :), ...
          interior_hist(boundary_hist_1(:,2)+1, :)];
        
        features = [features_boundary, features_interior];
%         if(size(features,2)~=linkage_model.n_dimension_insection)
%           error('Segment merge feature dimensions don''t match during training and testing.');
%         end
        weak_learner = tree_node_w(linkage_model.model_config.tree_depth_insection); %#ok<NASGU>
        confidence_boundary_1 = ...
          Classify(linkage_model.classifier_insection.GLearners, ...
          linkage_model.classifier_insection.GWeights, ...
          features');
        
        boundary_hist_1 = boundary_hist_1(...
          confidence_boundary_1>link_positive_threshold | ...
          confidence_boundary_1<0, :);
        confidence_boundary_1 = confidence_boundary_1(...
          confidence_boundary_1>link_positive_threshold | ...
          confidence_boundary_1<0);
        
        boundary_hist_2 = collect_segment_pair_stats_boundary_hist(...
          uint32(seg_2_original), uint8(255*images_2{tile_2}), ...
          uint8(linkage_model.feature_config.boundary_hist.bin_size));
        boundary_hist_2 = boundary_hist_2(min(boundary_hist_2(:,1:2), [], 2)>0, :);
        boundary_hist_2 = boundary_hist_2(...
          boundary_hist_2(:,1)~=boundary_hist_2(:,2), :);
        features_boundary = interval_sum(boundary_hist_2(:, 3:end));
        
        interior_hist_t = collect_segment_stats_interior_hist(...
          uint32(seg_2_original), uint8(255*images_2{tile_2}), ...
          uint8(feature_config.interior_hist.bin_size));
        interior_hist = [];
        interior_hist(interior_hist_t(:,1)+1, :) = interior_hist_t(:, 2:end); %#ok<AGROW>
        interior_hist = interval_sum(interior_hist);
        features_interior = [interior_hist(boundary_hist_2(:,1)+1, :), ...
          interior_hist(boundary_hist_2(:,2)+1, :)];
        
        features = [features_boundary, features_interior];
%         if(size(features,2)~=linkage_model.n_dimension_insection)
%           error('Segment merge feature dimensions don''t match during training and testing.');
%         end
        confidence_boundary_2 = ...
          Classify(linkage_model.classifier_insection.GLearners, ...
          linkage_model.classifier_insection.GWeights, ...
          features');
        
        boundary_hist_2 = boundary_hist_2(...
          confidence_boundary_2>link_positive_threshold | ...
          confidence_boundary_2<0, :);
        confidence_boundary_2 = confidence_boundary_2(...
          confidence_boundary_2>link_positive_threshold | ...
          confidence_boundary_2<0);
        
        min_node_id_1 = 1;
        max_node_id_1 = max(relabeling_1);
        min_node_id_2 = max_node_id_1 + 1;
        max_node_id_2 = max(relabeling_2);
        fprintf('min_node_id_1: %d, max_node_id_1: %d\n', ...
          min_node_id_1, max_node_id_1);
        fprintf('min_node_id_2: %d, max_node_id_2: %d\n', ...
          min_node_id_2, max_node_id_2);
        
        n_node = max_node_id_2;
        fprintf('n_node: %d\n', n_node);
        
        fprintf('building F matrix for co-clustering\n');
        F = eye(n_node, n_node);
        
        F(sub2ind(size(F), boundary_hist_1(:,1), boundary_hist_1(:,2))) = ...
          confidence_boundary_1;
        F(sub2ind(size(F), boundary_hist_1(:,2), boundary_hist_1(:,1))) = ...
          confidence_boundary_1;
        
        F(sub2ind(size(F), boundary_hist_2(:,1), boundary_hist_2(:,2))) = ...
          confidence_boundary_2;
        F(sub2ind(size(F), boundary_hist_2(:,2), boundary_hist_2(:,1))) = ...
          confidence_boundary_2;
        
        F(sub2ind(size(F), links_temp(:,1), links_temp(:,2))) = ...
          links_temp(:,3);
        F(sub2ind(size(F), links_temp(:,2), links_temp(:,1))) = ...
          links_temp(:,3);
        
        if(isfield(model_config, 'LP') && isfield(model_config.LP, 'all_tri_ineq') && ...
            model_config.LP.all_tri_ineq)
          A = ones(size(F));
        else
          A = double(F~=0);
        end
        
        lambda =  -1;
        fprintf('lambda: %g\n', lambda);
        Lambda = zeros(size(F))+lambda;
        
        G = F*100 + Lambda;
        G(G==0) = lambda;
        
        fprintf('Dumping linear program.\n');
        lp_file_name = [config_global.temp_dir, 'joint_link_lp3.lp'];
        dump_linear_program_p3_v2(double(G), double(A), lp_file_name);
        
        fprintf('Optimizing.\n');
        lp_sol_file_name = [config_global.temp_dir, 'joint_link_lp3.sol'];
        lp_solv_file_name = [config_global.temp_dir, 'joint_link_lp3.solv'];
        tic
        [status, result] = ...
          system(['lp_solve -lp ', lp_file_name, ' > ', lp_sol_file_name]);
        if(status~=0)
          error(result);
        end
        [status, result] = system(['cat ', lp_sol_file_name, ' | grep "x_" | ', ...
          'awk ''{print $0;}'' > ', lp_solv_file_name]);
        if(status~=0)
          error(result);
        end
        toc
        
        R = eye(n_node, n_node);
        fin = fopen(lp_solv_file_name, 'rt');
        while(true)
          v_name = fscanf(fin, '%s', [1 1]);
          if(feof(fin)==1)
            break;
          end
          v_value = 1-fscanf(fin, '%g', 1);
          if(v_value<=0)
            continue;
          end
          
          v_name_brkt = strfind(v_name, '_');
          node_i = str2double(v_name(v_name_brkt(1)+1:v_name_brkt(2)-1)) + 1;
          node_j = str2double(v_name(v_name_brkt(2)+1:end)) + 1;
          R(node_i, node_j) = v_value;
          R(node_j, node_i) = v_value;
        end
        fclose(fin);
        
        if(nnz((R>0) & (A==0))~=0)
          error('Spurious variables in LP result');
        end
        
        cc_id = get_connected_components_c(sparse(R), 0.95);
        links_temp = [...
          label_pairs_original(:,[1 2]), ...
          cc_id(label_pairs_original(:,1))==cc_id(label_pairs_original(:,2))];
      end
      
      if(~isempty(links_boost_area))
        links_temp = [links_temp; links_boost_area(links_boost_area(:,3)>2,:)]; %#ok<AGROW>
      end
      
      if(isfield(linkage_config, 'modified_segment') || ...
          isfield(apply_config, 'retain_links_subset'))
        links_temp = [links_temp; ...
          label_pairs_original, zeros(size(label_pairs_original,1),1)-2^14]; %#ok<AGROW>
      end

      if(~isempty(links_temp))
        links_temp = links_temp(links_temp(:,1)>0 & links_temp(:,2)>0, :);
        links_temp = [label_backref_1(1+links_temp(:,1))', ...
          label_backref_2(1+links_temp(:,2))', links_temp(:,3)];
      end
      
      if(~isempty(links_temp))
        links_temp = unique(links_temp, 'rows');
        
        % remove locked segments from the linkage graph
        if(isfield(segment_2D_label_map_1{tile_1}, 'locked_labels'))
          to_retain = ~ismember(links_temp(:,1), ...
            segment_2D_label_map_1{tile_1}.locked_labels);
          links_temp = links_temp(to_retain, :);
        end
        if(isfield(segment_2D_label_map_2{tile_2}, 'locked_labels'))
          to_retain = ~ismember(links_temp(:,2), ...
            segment_2D_label_map_2{tile_2}.locked_labels);
          links_temp = links_temp(to_retain, :);
        end
        
        links_temp = links_temp(links_temp(:,1)~=0 & links_temp(:,2)~=0, :);
      end
      
      links_3D_p = links_temp;
      fprintf('number of overlap pairs: %d\n', size(links_3D_p,1));
      
      % Include proofread linkage graph if asked to do so.
      if(isfield(linkage_config, 'include_proofread') && ...
          linkage_config.include_proofread.is_enabled)
        fprintf('Including proofread linkage graph, proofread_name: %s\n', ...
          linkage_config.include_proofread.proofread_name);
        links_3D_p = include_proofread_linkage_graph(...
          image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, ...
          linkage_config, links_3D_p);
      end
      
      if(~isempty(links_3D_p))
        fprintf('saving linkage graph:\n%s\n', file_name_linkage_graph);
        save2(file_name_linkage_graph, 'links_3D_p');
      end
    end
  end

  images_1 = images_2;  
  image_prefixes_1 = image_prefixes_2;
  segment_2D_label_map_1 = segment_2D_label_map_2;
  is_to_be_processed_1 = is_to_be_processed_2;
  fprintf('\n');
end;

if(linkage_config.is_verbose)
  fprintf('STOP: gen_link_gph_overlap_hist_boundary_intr_hist_lp_wekarf_dmesh\n');
end
return
end
