function gen_linkage_gph_overlap_hist_boost_deformable_mesh(config)
% gen_linkage_gph_intensity_pair_v2_boost_deformable_mesh(config)
% generate linkage graph using boosted classifier trained on
% intensity histograms features
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
% v0  ~03202008 init code
% v1  04112008  modified for reconstruction pipeline
% v2  07152008  appended multiresolution pyramid of the 2D histograms.
% v3  12122008  modified for deformable mesh
%

global config_global

stack_config = config.stack;
linkage_config = config.linkage;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
apply_config = linkage_config.apply;
dmesh_config = config.align.linkage_align.deformable_mesh;

version_names{1} = 'overlap_hist';
version_names{2} = 'overlap_hist_shrink_LS';
version_id = find(strcmp(feature_config.type, version_names)); %#ok<EFIND>
if(isempty(version_id))
  error('Feature type does not match with called function. Exiting');
end;
if(strcmp(model_config.type, 'boost')==0)
  error('Classifier type does not match with called function. Exiting');
end;

get_linkage_dir(config);

dmesh_dir = get_deformable_mesh_dir(config);

junk = tree_node_w(1); %#ok<NASGU>
linkage_model = load2([get_reconstruction_dir(config), linkage_config.dir, 'link_model', '.', ...
  model_config.type, model_config.suffix, '.', ...
  feature_config.type, feature_config.suffix, ...
  apply_config.model_suffix, '.mat'], ...
  'classifier', 'n_dimension', 'model_config', 'feature_config');
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
      
      fprintf('Computing linkage confidences\n');
      if(isfield(linkage_config, 'modified_segment') || ...
          isfield(apply_config, 'retain_links_subset'))
        [features, label_pairs, label_pairs_original] = ...
          get_linkage_features_overlap_hist_shrink_LS(...
          segment_2D_label_map_1{tile_1}.label_map, images_1{tile_1}, ...
          segment_2D_label_map_2{tile_2}.label_map, images_2{tile_2}, ...
          transforms_tp, transforms_tp_rev, feature_config, ...
          segment_2D_label_map_1{tile_1}.original_label_map, ...
          segment_2D_label_map_2{tile_2}.original_label_map);
        if(~isempty(label_pairs_original))
          label_pairs_original = label_pairs_original(...
            label_pairs_original(:,1)>0 & label_pairs_original(:,2)>0, :);
        end
      else
        [features, label_pairs] = get_linkage_features_overlap_hist_shrink_LS(...
          segment_2D_label_map_1{tile_1}.label_map, images_1{tile_1}, ...
          segment_2D_label_map_2{tile_2}.label_map, images_2{tile_2}, ...
          transforms_tp, transforms_tp_rev, feature_config);
      end
      if(~isempty(features))
        link_confidence = ...
          Classify(linkage_model.classifier.GLearners, linkage_model.classifier.GWeights, ...
          features');
        links_temp = [label_pairs, link_confidence'];
      end

      if(~isempty(links_temp))
        links_temp = sortrows(links_temp);
        
        if(isfield(apply_config, 'retain_links_subset'))
          switch(apply_config.retain_links_subset)
            case 'max_confidence'
              l1 = links_temp;
              l1 = sortrows(l1, 3);
              [junk, I] = unique(l1(:,1), 'last'); %#ok<ASGLU>
              l1 = l1(I,:);
              
              l2 = links_temp;
              l2 = sortrows(l2, 3);
              [junk, I] = unique(l2(:,2), 'last'); %#ok<ASGLU>
              l2 = l2(I,:);
              
              l3 = links_temp(links_temp(:,3)>2.5, :);
              
              links_temp = [l1; l2; l3];
              links_temp = unique(links_temp, 'rows');
          end
        end
      end
      
      [features, label_pairs] = get_linkage_features_overlap_area_shrink_LS(...
        segment_2D_label_map_1{tile_1}.label_map, images_1{tile_1}, ...
        segment_2D_label_map_2{tile_2}.label_map, images_2{tile_2}, ...
        transforms_tp, transforms_tp_rev);
      if(~isempty(features))
        link_confidence = (features(:,1)>10000).*(1 + features(:,1)/10000);
        links_temp = [links_temp; label_pairs, link_confidence]; %#ok<AGROW>
      end
      
      if(false && ~( isempty(segment_2D_label_map_1{tile_1}.label_map) || ...
          isempty(segment_2D_label_map_2{tile_2}.label_map) || ...
          isempty(transforms_tp) || ...
          isempty(transforms_tp.transforms) || ...
          isempty(transforms_tp.map_mask) || ...
          size(transforms_tp.transforms, 2)<=4 || ...
          size(transforms_tp_rev.transforms, 2)<=4))
        [coords_x_1, coords_y_1] = meshgrid(1:0.5:size(images_1{tile_1},2), ...
          1:0.5:size(images_1{tile_1},1));
        coords_x_1 = coords_x_1(:)';
        coords_y_1 = coords_y_1(:)';
        
        tm = transforms_tp.map_mask;
        tm(tm>0) = tm(tm>0) - config_global.TRANSFORMATION_ID_OFFSET + 1;
        t1 = [0 transforms_tp.transforms(1,:)];
        t2 = [0 transforms_tp.transforms(2,:)];
        t3 = [0 transforms_tp.transforms(3,:)];
        t4 = [0 transforms_tp.transforms(4,:)];
        t5 = [0 transforms_tp.transforms(5,:)];
        t6 = [0 transforms_tp.transforms(6,:)];
        tid = tm(sub2ind(size(tm), round(coords_y_1), round(coords_x_1)))'+1;
        coords_x_1_t = t1(tid).*coords_x_1 + t3(tid).*coords_y_1 + t5(tid);
        coords_y_1_t = t2(tid).*coords_x_1 + t4(tid).*coords_y_1 + t6(tid);

        coords_x_1 = round(coords_x_1);
        coords_y_1 = round(coords_y_1);
        coords_x_1_t = round(coords_x_1_t);
        coords_y_1_t = round(coords_y_1_t);
        
        to_retain = coords_x_1_t>=1 & coords_x_1_t<=size(images_2{tile_2},2) &...
          coords_y_1_t>=1 & coords_y_1_t<=size(images_2{tile_2},1);
        coords_x_1 = coords_x_1(to_retain);
        coords_y_1 = coords_y_1(to_retain);
        coords_x_1_t = coords_x_1_t(to_retain);
        coords_y_1_t = coords_y_1_t(to_retain);
        
        if(false && isfield(segment_2D_label_map_1{tile_1}, 'original_label_map'))
          segment_0 = zeros(size(images_2{tile_2}));
          segment_0(sub2ind(size(images_2{tile_2}), coords_y_1_t, coords_x_1_t)) = ...
            segment_2D_label_map_1{tile_1}.original_label_map(sub2ind(size(images_1{tile_1}), ...
            coords_y_1, coords_x_1));
          segment_1 = segment_2D_label_map_2{tile_2}.original_label_map;
        else
          segment_0 = zeros(size(images_2{tile_2}));
          segment_0(sub2ind(size(images_2{tile_2}), coords_y_1_t, coords_x_1_t)) = ...
            segment_2D_label_map_1{tile_1}.label_map(sub2ind(size(images_1{tile_1}), ...
            coords_y_1, coords_x_1));
          segment_1 = segment_2D_label_map_2{tile_2}.label_map;
        end
        
        border_mask_0 = imdilate(segment_0>0, strel('disk', 2));
        border_mask_1 = imdilate(segment_1>0, strel('disk', 2));
        stitch_segments_0 = unique(segment_0(:));
        
        c.align_segmentation.is_verbose_figures = false;
        c.align_segmentation.delta_margin = 7;
        c.align_segmentation.min_n_votes = 40;

        [junk, segment_correspondence] = ...
          correspond_segment_boundary_bipartite_match...
          (segment_0, border_mask_0, stitch_segments_0, segment_1, ...
          border_mask_1, c); %#ok<ASGLU>
        segment_correspondence = segment_correspondence';
        segment_correspondence = segment_correspondence(...
          segment_correspondence(:,1)>0 & segment_correspondence(:,2)>0, :);
        l1 = segment_correspondence;
        l1 = sortrows(l1, 3);
        [junk, I] = unique(l1(:,1), 'last'); %#ok<ASGLU>
        l1 = l1(I,:);
        
        l2 = segment_correspondence;
        l2 = sortrows(l2, 3);
        [junk, I] = unique(l2(:,2), 'last'); %#ok<ASGLU>
        l2 = l2(I,:);
        
        l12 = [l1; l2];
        
        links_temp = [links_temp; ...
          l12(:,[1 2]), 3.1416*(l12(:,3)>c.align_segmentation.min_n_votes)]; %#ok<AGROW>
      end
      
      if(isfield(linkage_config, 'modified_segment') || ...
          isfield(apply_config, 'retain_links_subset'))
        links_temp = [links_temp; ...
          label_pairs_original, zeros(size(label_pairs_original,1),1)-2^14]; %#ok<AGROW>
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

fprintf('done\n');
