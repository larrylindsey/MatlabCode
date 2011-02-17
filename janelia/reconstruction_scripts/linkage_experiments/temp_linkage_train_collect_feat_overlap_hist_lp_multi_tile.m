function linkage_train_collect_feat_overlap_hist_lp_multi_tile(config)
% linkage_train_collect_feat_overlap_hist_multi_tile(config)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v1  07092008  init code
% v2  06222009  split data section-wise
% v3  07282009  include modified segmentations.
%

if(~isfield(config.linkage, 'is_verbose'))
  config.linkage.is_verbose = true;
end
if(~isfield(config.linkage, 'is_verbose_figures'))
  config.linkage.is_verbose_figures = false;
end
if(config.linkage.is_verbose)
  fprintf('START: linkage_train_collect_feat_overlap_hist_lp_multi_tile\n');
end

stack_config = config.stack;
linkage_config = config.linkage;
feature_config = linkage_config.feature;
dmesh_config = config.align.linkage_align.deformable_mesh;

version_names{1} = 'overlap_hist_boundary_interior_hist';
version_names{2} = 'overlap_hist_boundary_interior_hist_shrink_LS';
version_id = find(strcmp(feature_config.type, version_names)); %#ok<EFIND>
if(isempty(version_id))
  error('Feature type does not match with called function. Exiting');
end;

get_linkage_dir(config);
dmesh_dir = get_deformable_mesh_dir(config);
seg_dir_gt = [get_reconstruction_dir(config), config.segmentation_2D(end).dir, ...
  linkage_config.train.segment_method_ground_truth, '/'];

if(~isfield(feature_config, 'suffix'))
  feature_config.suffix = '';
end

map_param.normalized_area_threshold = 0.8;

case_id = stack_config.case_ids(1);
fprintf('case_id: %d\n', case_id);
[images_1, image_prefixes_1, image_sub_dirs_1, is_to_be_processed_1, size_images_1] = ...
  get_image_from_stack(config, case_id); %#ok<ASGLU>
if(isfield(feature_config, 'filter_version') && ...
    ~isempty(feature_config.filter_version))
  for tile = 1:length(images_1)
    images_1{tile} = filter_image2(images_1{tile}, feature_config.filter_version, config);
  end
end
segment_2D_label_map_1 = load_segment_maps(image_prefixes_1, size_images_1, config);
if(linkage_config.is_verbose)
  fprintf('Loading ground-truth segment maps\n');
end
seg_gt_1 = {};
for tile = 1:length(image_prefixes_1)
  seg_gt_1{tile} = load2([seg_dir_gt, image_prefixes_1{tile}, ...
    linkage_config.train.segment_suffix_ground_truth,'.mat']); %#ok<AGROW>
  if(~isfield(seg_gt_1{tile}, 'label_map'))
    sp_dir_gt = [get_reconstruction_dir(config), config.segmentation_2D(end).dir, ...
      linkage_config.train.superpixel_method_ground_truth, '/'];
    sp_gt_1 = load2([sp_dir_gt, image_prefixes_1{tile}, ...
      linkage_config.train.superpixel_suffix_ground_truth,'.mat']);
    seg_gt_1{tile}.label_map = apply_mapping(sp_gt_1.label_map, ...
      seg_gt_1{tile}.superpixel_to_seg_label); %#ok<AGROW>
  end
end

for i = 2:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  fprintf('case_id: %d\n', case_id);

  [images_2, image_prefixes_2, image_sub_dirs_2, ...
    is_to_be_processed_2, size_images_2] = get_image_from_stack(config, case_id); %#ok<ASGLU>
  
  if(isfield(feature_config, 'filter_version') && ...
      ~isempty(feature_config.filter_version))
    for tile = 1:length(image_prefixes_2)
      images_2{tile} = filter_image2(images_2{tile}, feature_config.filter_version, config);
    end
  end
  segment_2D_label_map_2 = load_segment_maps(image_prefixes_2, size_images_2, config);
  if(linkage_config.is_verbose)
    fprintf('Loading ground-truth segment maps\n');
  end
  seg_gt_2 = {};
  for tile = 1:length(image_prefixes_2)
    seg_gt_2{tile} = load2([seg_dir_gt, image_prefixes_2{tile}, ...
      linkage_config.train.segment_suffix_ground_truth,'.mat']); %#ok<AGROW>
    if(~isfield(seg_gt_2{tile}, 'label_map'))
      sp_dir_gt = [get_reconstruction_dir(config), config.segmentation_2D(end).dir, ...
        linkage_config.train.superpixel_method_ground_truth, '/'];
      sp_gt_2 = load2([sp_dir_gt, image_prefixes_2{tile}, ...
        linkage_config.train.superpixel_suffix_ground_truth,'.mat']);
      seg_gt_2{tile}.label_map = apply_mapping(sp_gt_2.label_map, ...
        seg_gt_2{tile}.superpixel_to_seg_label); %#ok<AGROW>
    end
  end

  
  % collect linkage features of pairs of overlapping segments for each pair
  % of overlapping tiles
  for tile_1 = 1:length(image_prefixes_1)
    for tile_2 = 1:length(image_prefixes_2)
      fprintf('--- Tile pair ---\n%d,%d\n', tile_1, tile_2);
      fprintf('%s\n%s\n', image_prefixes_1{tile_1}, image_prefixes_2{tile_2});
      
      if(is_to_be_processed_1(tile_1)==0 && is_to_be_processed_2(tile_2)==0)
        fprintf('Both tiles have is_to_be_processed=false, skipping this pair\n');
        continue;
      end
      
      file_name_prefix = get_file_name_from_tuple(...
        [get_reconstruction_dir(config), linkage_config.dir], ...
        image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'flp.');
      file_name_suffix = ['.', feature_config.type, feature_config.suffix];
      file_name_features = [file_name_prefix, file_name_suffix, ...
        config.segmentation_choose.choice.seg_suffix, '.mat'];
      fprintf('Deleting the target feature file if it already exists\n');
      delete(get_storage_file_name(file_name_features));
      
      if(isempty(seg_gt_1{tile_1}.label_map) || ...
          max(seg_gt_1{tile_1}.label_map(:))<=0)
        if(linkage_config.is_verbose)
          fprintf('Ground-truth segment map of tile 1 is empty, skipping\n');
        end
        continue;
      end
      if(isempty(seg_gt_2{tile_2}.label_map) || ...
          max(seg_gt_2{tile_2}.label_map(:))<=0)
        if(linkage_config.is_verbose)
          fprintf('Ground-truth segment map of tile 2 is empty, skipping\n');
        end
        continue;
      end
      
      file_name_suffix = get_file_name_from_tuple(get_linkage_dir(config), ...
        image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'lkp.');
      file_name = [file_name_suffix, ...
        linkage_config.train.linkage_suffix_ground_truth, '.mat'];
      try
        if(linkage_config.is_verbose)
          fprintf('Attempting to load ground-truth:\n%s\n', file_name);
        end
        linkage_graph_gt = load2(file_name, 'links_3D_p');
      catch %#ok<CTCH>
        if(linkage_config.is_verbose)
          fprintf('failed, skipping.\n');
        end
        continue;
      end
      if(linkage_config.is_verbose)
        fprintf('loaded.\n');
      end
      
      linkage_graph_gt.links_3D_p = linkage_graph_gt.links_3D_p(...
        linkage_graph_gt.links_3D_p(:,1)>0 & linkage_graph_gt.links_3D_p(:,2)>0, :);
      if(isempty(linkage_graph_gt.links_3D_p))
        if(linkage_config.is_verbose)
          fprintf('Ground-truth linkage graph empty, skipping.\n');
        end
        continue;
      end
      
      fprintf('Loading transforms between tile 1 and tile 2\n');
      [transforms_tp, transforms_tp_rev] = ...
        load_tile_pair_deformable_mesh_transforms(image_prefixes_1{tile_1}, ...
        image_prefixes_2{tile_2}, dmesh_config, dmesh_dir);
      
      fprintf('Computing linkage features\n');
      if(isfield(linkage_config, 'modified_segment'))
        [features, seg_label_pairs] = ...
          get_linkage_features_overlap_hist_shrink_LS(...
          segment_2D_label_map_1{tile_1}.label_map, images_1{tile_1}, ...
          segment_2D_label_map_2{tile_2}.label_map, images_2{tile_2}, ...
          transforms_tp, transforms_tp_rev, feature_config, ...
          segment_2D_label_map_1{tile_1}.original_label_map, ...
          segment_2D_label_map_2{tile_2}.original_label_map);
      else
        [features, seg_label_pairs] = get_linkage_features_overlap_hist_shrink_LS(...
          segment_2D_label_map_1{tile_1}.original_label_map, images_1{tile_1}, ...
          segment_2D_label_map_2{tile_2}.original_label_map, images_2{tile_2}, ...
          transforms_tp, transforms_tp_rev, feature_config);
      end
      if(isempty(features))
        fprintf('No feature vectors, skipping.\n');
        continue;
      end
      features = features(seg_label_pairs(:,1)>0 & seg_label_pairs(:,2)>0, :);
      seg_label_pairs = ...
        seg_label_pairs(seg_label_pairs(:,1)>0 & seg_label_pairs(:,2)>0, :);
      if(isempty(features))
        fprintf('No feature vectors, skipping.\n');
        continue;
      end
      
      if(linkage_config.is_verbose)
        fprintf('Mapping automatic segments to ground-truth bodies\n');
      end
      segment_offset = max(linkage_graph_gt.links_3D_p(:,1)) + 1;
      linkage_graph_gt.links_3D_p(:,2) = ...
        linkage_graph_gt.links_3D_p(:,2) + segment_offset;
      seg_gt_2_temp = seg_gt_2{tile_2};
      seg_gt_2_temp.label_map(seg_gt_2_temp.label_map>0) = ...
        seg_gt_2_temp.label_map(seg_gt_2_temp.label_map>0) + segment_offset;
      sec_seg_2_body_map = get_sec_seg_2_body_map_from_links_3D_with_dummy(...
        {linkage_graph_gt.links_3D_p}, 0.5);
      seg_2_body_map_temp = sec_seg_2_body_map(:, 2:3);
      body_map_gt_1 = apply_mapping(seg_gt_1{tile_1}.label_map, ...
        seg_2_body_map_temp);
      body_map_gt_2 = apply_mapping(seg_gt_2_temp.label_map, ...
        seg_2_body_map_temp);
      
      seg_map_auto_2_gt_1 = assign_auto_seg_label_to_groundtruth(...
        segment_2D_label_map_1{tile_1}.original_label_map, body_map_gt_1, map_param);
      seg_map_auto_2_gt_2 = assign_auto_seg_label_to_groundtruth(...
        segment_2D_label_map_2{tile_2}.original_label_map, body_map_gt_2, map_param);
      
      to_retain = apply_mapping(seg_label_pairs(:,1), seg_map_auto_2_gt_1)>0 & ...
        apply_mapping(seg_label_pairs(:,2), seg_map_auto_2_gt_2)>0;
      features = features(to_retain, :); %#ok<NASGU>
      seg_label_pairs = seg_label_pairs(to_retain, :);
      
      link_labels = 2*(apply_mapping(seg_label_pairs(:,1), seg_map_auto_2_gt_1)...
        ==apply_mapping(seg_label_pairs(:,2), seg_map_auto_2_gt_2))-1; %#ok<NASGU>
      
      link_pair_ids = seg_label_pairs; %#ok<NASGU>
      
      if(linkage_config.is_verbose)
        fprintf('Collecting within-section link features\n');
      end
      boundary_hist_1 = collect_segment_pair_stats_boundary_hist(...
        uint32(segment_2D_label_map_1{tile_1}.original_label_map), ...
        uint8(255*images_1{tile_1}), ...
        uint8(feature_config.boundary_hist.bin_size));
      boundary_hist_1 = boundary_hist_1(boundary_hist_1(:,1)>0 & ...
        boundary_hist_1(:,2)>0, :);
      boundary_hist_1 = boundary_hist_1(boundary_hist_1(:,1)~=...
        boundary_hist_1(:,2), :);
      boundary_hist_1 = boundary_hist_1(...
        apply_mapping(boundary_hist_1(:,1), seg_map_auto_2_gt_1)>0 & ...
        apply_mapping(boundary_hist_1(:,2), seg_map_auto_2_gt_1)>0, :);
      
      boundary_feature_1 = interval_sum(boundary_hist_1(:, 3:end));
      
      interior_hist_t = collect_segment_stats_interior_hist(...
        uint32(segment_2D_label_map_1{tile_1}.original_label_map), ...
        uint8(255*images_1{tile_1}), ...
        uint8(feature_config.interior_hist.bin_size));
      interior_hist_t = interior_hist_t(interior_hist_t(:,1)>0, :);
      interior_hist_1 = [];
      interior_hist_1(interior_hist_t(:,1)+1, :) = interior_hist_t(:, 2:end); %#ok<AGROW>
      interior_hist_1 = interval_sum(interior_hist_1);
      interior_feature_1 = [interior_hist_1(boundary_hist_1(:,1)+1, :), ...
        interior_hist_1(boundary_hist_1(:,2)+1, :); ...
        interior_hist_1(boundary_hist_1(:,2)+1, :), ...
        interior_hist_1(boundary_hist_1(:,1)+1, :)];

      boundary_hist_2 = collect_segment_pair_stats_boundary_hist(...
        uint32(segment_2D_label_map_2{tile_2}.original_label_map), ...
        uint8(255*images_2{tile_2}), ...
        uint8(feature_config.boundary_hist.bin_size));
      boundary_hist_2 = boundary_hist_2(boundary_hist_2(:,1)>0 & ...
        boundary_hist_2(:,2)>0, :);
      boundary_hist_2 = boundary_hist_2(boundary_hist_2(:,1)~=...
        boundary_hist_2(:,2), :);
      boundary_hist_2 = boundary_hist_2(...
        apply_mapping(boundary_hist_2(:,1), seg_map_auto_2_gt_2)>0 & ...
        apply_mapping(boundary_hist_2(:,2), seg_map_auto_2_gt_2)>0, :);
      boundary_feature_2 = interval_sum(boundary_hist_2(:, 3:end));
      
      interior_hist_t = collect_segment_stats_interior_hist(...
        uint32(segment_2D_label_map_2{tile_2}.original_label_map), ...
        uint8(255*images_2{tile_2}), ...
        uint8(feature_config.interior_hist.bin_size));
      interior_hist_t = interior_hist_t(interior_hist_t(:,1)>0, :);
      interior_hist_2 = [];
      interior_hist_2(interior_hist_t(:,1)+1, :) = interior_hist_t(:, 2:end); %#ok<AGROW>
      interior_hist_2 = interval_sum(interior_hist_2);
      interior_feature_2 = [interior_hist_2(boundary_hist_2(:,1)+1, :), ...
        interior_hist_2(boundary_hist_2(:,2)+1, :); ...
        interior_hist_2(boundary_hist_2(:,2)+1, :), ...
        interior_hist_2(boundary_hist_2(:,1)+1, :)];
      
      features_insection = [ ...
        [[boundary_feature_1; boundary_feature_1], interior_feature_1]; ...
        [[boundary_feature_2; boundary_feature_2], interior_feature_2]; ...
        ]; %#ok<NASGU>
      
      labels_insection_1 = repmat(...
        2*(apply_mapping(boundary_hist_1(:,1), seg_map_auto_2_gt_1)...
        ==apply_mapping(boundary_hist_1(:,2), seg_map_auto_2_gt_1))-1, [2 1]);
      labels_insection_2 = repmat(...
        2*(apply_mapping(boundary_hist_2(:,1), seg_map_auto_2_gt_2)...
        ==apply_mapping(boundary_hist_2(:,2), seg_map_auto_2_gt_2))-1, [2 1]);
      labels_insection = [labels_insection_1; labels_insection_2]; %#ok<NASGU>
      
      pair_ids_insection = [...
        repmat(boundary_hist_1(:,1:2), [2 1]); ...
        repmat(boundary_hist_2(:,1:2), [2 1])]; %#ok<NASGU>
      
      fprintf('saving linkage features:\n%s\n', file_name_features);
      save2(file_name_features, 'features', 'link_labels', ...
        'features_insection', 'labels_insection', ...
        'link_pair_ids', 'pair_ids_insection');
    end
  end

  images_1 = images_2;  
  image_prefixes_1 = image_prefixes_2;
  segment_2D_label_map_1 = segment_2D_label_map_2;
  is_to_be_processed_1 = is_to_be_processed_2;
  seg_gt_1 = seg_gt_2;
end;

if(linkage_config.is_verbose)
  fprintf('STOP: linkage_train_collect_feat_overlap_hist_lp_multi_tile\n');
end
return
end
