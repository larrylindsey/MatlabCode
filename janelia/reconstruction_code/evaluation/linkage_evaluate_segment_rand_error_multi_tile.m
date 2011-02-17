function linkage_evaluate_segment_rand_error_multi_tile(config)
% linkage_evaluate_segment_rand_error_multi_tile(config)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v1  07092008  init code
% v2  06222009  split data section-wise
% v3  07282009  include modified segmentations.
%

if(~isfield(config.linkage.evaluate, 'is_verbose'))
  config.linkage.evaluate.is_verbose = true;
end
if(~isfield(config.linkage.evaluate, 'is_verbose_figures'))
  config.linkage.evaluate.is_verbose_figures = false;
end
if(config.linkage.evaluate.is_verbose)
  fprintf('START: linkage_evaluate_segment_rand_error_multi_tile\n');
end

stack_config = config.stack;
linkage_config = config.linkage;
evaluate_config = linkage_config.evaluate;

get_linkage_dir(config);
evaluate_save_dir = [get_linkage_dir(config), config.linkage.evaluate.save_dir];
check_for_dir(evaluate_save_dir);

image_prefixes_2 = {};
map_param.normalized_area_threshold = 0.8;
for layer_id = 1:length(stack_config.case_ids)-1
  case_id = stack_config.case_ids(layer_id);
  case_id_1 = stack_config.case_ids(layer_id+1);
  fprintf('case_id: %d, case_id_1: %d\n', case_id, case_id_1);
  
  if(isempty(image_prefixes_2))
    image_prefixes_1 = get_image_prefixes_subdirs(config, case_id);
  else
    image_prefixes_1 = image_prefixes_2;
  end
  image_prefixes_2 = get_image_prefixes_subdirs(config, case_id_1);
  
  seg_dir_gt = [get_reconstruction_dir(config), config.segmentation_2D(end).dir, ...
    evaluate_config.segment_method_ground_truth, '/'];
  for tile_1 = 1:length(image_prefixes_1)
    for tile_2 = 1:length(image_prefixes_2)
      if(evaluate_config.is_verbose)
        fprintf('tile 1:%d, tile 2:%d\n', tile_1, tile_2);
        fprintf('tile 1: %s\ntile 2: %s\n', image_prefixes_1{tile_1}, ...
          image_prefixes_2{tile_2});
      end
      
      file_name_prefix = get_file_name_from_tuple(evaluate_save_dir, ...
        image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'elp.');
      file_name_suffix = [evaluate_config.linkage_suffix_automatic, ...
        evaluate_config.linkage_suffix_ground_truth];
      file_name_linkage_eval = [file_name_prefix, file_name_suffix, '.mat'];
      fprintf('Deleting the linkage evaluation file if it already exists\n');
      delete(get_storage_file_name(file_name_linkage_eval));

      file_name_suffix = get_file_name_from_tuple(get_linkage_dir(config), ...
        image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'lkp.');
      file_name = [file_name_suffix, ...
        evaluate_config.linkage_suffix_automatic, '.mat'];
      try
        if(evaluate_config.is_verbose)
          fprintf('Attempting to load automatic linkage:\n%s\n', file_name);
        end
        linkage_graph_auto = load2(file_name, 'links_3D_p');
      catch %#ok<CTCH>
        if(evaluate_config.is_verbose)
          fprintf('failed, skipping.\n');
        end
        continue;
      end
      if(evaluate_config.is_verbose)
        fprintf('loaded.\n');
      end
      if(isempty(linkage_graph_auto.links_3D_p))
        if(evaluate_config.is_verbose)
          fprintf('Automatic linkage graph empty, skipping.\n');
        end
        continue;
      end
      
      file_name_suffix = get_file_name_from_tuple(get_linkage_dir(config), ...
        image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'lkp.');
      file_name = [file_name_suffix, ...
        evaluate_config.linkage_suffix_ground_truth, '.mat'];
      try
        if(evaluate_config.is_verbose)
          fprintf('Attempting to load ground-truth:\n%s\n', file_name);
        end
        linkage_graph_gt = load2(file_name, 'links_3D_p');
      catch %#ok<CTCH>
        if(evaluate_config.is_verbose)
          fprintf('failed, skipping.\n');
        end
        continue;
      end
      if(evaluate_config.is_verbose)
        fprintf('loaded.\n');
      end
      
      linkage_graph_gt.links_3D_p = linkage_graph_gt.links_3D_p(...
        linkage_graph_gt.links_3D_p(:,1)>0 & linkage_graph_gt.links_3D_p(:,2)>0, :);
      if(isempty(linkage_graph_gt.links_3D_p))
        if(evaluate_config.is_verbose)
          fprintf('Ground-truth linkage graph empty, skipping.\n');
        end
        continue;
      end
      
      if(evaluate_config.is_verbose)
        fprintf('Loading ground-truth segment for tile 1\n');
      end
      seg_gt_1 = load2([seg_dir_gt, image_prefixes_1{tile_1}, ...
        evaluate_config.segment_suffix_ground_truth,'.mat']);
      if(~isfield(seg_gt_1, 'label_map'))
        sp_dir_gt = [get_reconstruction_dir(config), config.segmentation_2D(end).dir, ...
          evaluate_config.superpixel_method_ground_truth, '/'];
        sp_gt_1 = load2([sp_dir_gt, image_prefixes_1{tile_1}, ...
          evaluate_config.superpixel_suffix_ground_truth,'.mat']);
        seg_gt_1.label_map = apply_mapping(sp_gt_1.label_map, ...
          seg_gt_1.superpixel_to_seg_label);
      end
      if(isempty(seg_gt_1.label_map) || max(seg_gt_1.label_map(:))<=0)
        if(evaluate_config.is_verbose)
          fprintf('Ground-truth segment map is empty, skipping\n');
        end
        continue;
      end
      
      if(evaluate_config.is_verbose)
        fprintf('Loading automatic segment for tile 1\n');
      end
      [seg_method, seg_suffix] = get_segmentation_suffix(config, image_prefixes_1{tile_1});
      seg_dir_auto = [get_reconstruction_dir(config), config.segmentation_2D(end).dir, ...
        seg_method, '/'];
      seg_auto_1 = load2([seg_dir_auto, image_prefixes_1{tile_1}, ...
        seg_suffix,'.mat']);
      if(~isfield(seg_auto_1, 'label_map'))
        if(isfield(config, 'segmentation_2D'))
          config_segmentation_2D = config.segmentation_2D;
          replace_flag = true;
        end
        config.segmentation_2D = config.superpixel_2_seg(1);
        [sp_method, sp_suffix] = ...
          get_superpixel_suffixes(config, image_prefixes_1{tile_1});
        if(replace_flag)
          config.segmentation_2D = config_segmentation_2D;
        end
        sp_dir_auto = [get_reconstruction_dir(config), config.segmentation_2D(end).dir, ...
          sp_method, '/'];
        sp_auto_1 = load2([sp_dir_auto, image_prefixes_1{tile_1}, ...
          sp_suffix,'.mat']);
        seg_auto_1.label_map = apply_mapping(sp_auto_1.label_map, ...
          seg_auto_1.superpixel_to_seg_label);
      end
      if(isempty(seg_auto_1.label_map) || max(seg_auto_1.label_map(:))<=0)
        fprintf('Automatic segment map is empty, skipping\n');
        continue;
      end
      
      if(evaluate_config.is_verbose)
        fprintf('Loading ground-truth segment for tile 2\n');
      end
      seg_gt_2 = load2([seg_dir_gt, image_prefixes_2{tile_2}, ...
        evaluate_config.segment_suffix_ground_truth,'.mat']);
      if(~isfield(seg_gt_2, 'label_map'))
        sp_dir_gt = [get_reconstruction_dir(config), config.segmentation_2D(end).dir, ...
          evaluate_config.superpixel_method_ground_truth, '/'];
        sp_gt_2 = load2([sp_dir_gt, image_prefixes_2{tile_2}, ...
          evaluate_config.superpixel_suffix_ground_truth,'.mat']);
        seg_gt_2.label_map = apply_mapping(sp_gt_2.label_map, ...
          seg_gt_2.superpixel_to_seg_label);
      end
      if(isempty(seg_gt_2.label_map) || max(seg_gt_2.label_map(:))<=0)
        if(evaluate_config.is_verbose)
          fprintf('Ground-truth segment map is empty, skipping\n');
        end
        continue;
      end
      
      if(evaluate_config.is_verbose)
        fprintf('Loading automatic segment for tile 2\n');
      end
      [seg_method, seg_suffix] = get_segmentation_suffix(config, image_prefixes_2{tile_2});
      seg_dir_auto = [get_reconstruction_dir(config), config.segmentation_2D(end).dir, ...
        seg_method, '/'];
      seg_auto_2 = load2([seg_dir_auto, image_prefixes_2{tile_2}, ...
        seg_suffix,'.mat']);
      if(~isfield(seg_auto_2, 'label_map'))
        if(isfield(config, 'segmentation_2D'))
          config_segmentation_2D = config.segmentation_2D;
          replace_flag = true;
        end
        config.segmentation_2D = config.superpixel_2_seg(1);
        [sp_method, sp_suffix] = ...
          get_superpixel_suffixes(config, image_prefixes_2{tile_2});
        if(replace_flag)
          config.segmentation_2D = config_segmentation_2D;
        end
        sp_dir_auto = [get_reconstruction_dir(config), config.segmentation_2D(end).dir, ...
          sp_method, '/'];
        sp_auto_2 = load2([sp_dir_auto, image_prefixes_2{tile_2}, ...
          sp_suffix,'.mat']);
        seg_auto_2.label_map = apply_mapping(sp_auto_2.label_map, ...
          seg_auto_2.superpixel_to_seg_label);
      end
      if(isempty(seg_auto_2.label_map) || max(seg_auto_2.label_map(:))<=0)
        fprintf('Automatic segment map is empty, skipping\n');
        continue;
      end
      
      if(evaluate_config.is_verbose)
        fprintf('Mapping automatic segments to ground-truth bodies\n');
      end
      segment_offset = max(linkage_graph_gt.links_3D_p(:,1)) + 1;
      linkage_graph_gt.links_3D_p(:,2) = ...
        linkage_graph_gt.links_3D_p(:,2) + segment_offset;
      seg_gt_2.label_map(seg_gt_2.label_map>0) = ...
        seg_gt_2.label_map(seg_gt_2.label_map>0) + segment_offset;
      sec_seg_2_body_map = get_sec_seg_2_body_map_from_links_3D_with_dummy(...
        {linkage_graph_gt.links_3D_p}, 0.5);
      seg_2_body_map_temp = sec_seg_2_body_map(:, 2:3);
      body_map_gt_1 = apply_mapping(seg_gt_1.label_map, ...
        seg_2_body_map_temp);
      body_map_gt_2 = apply_mapping(seg_gt_2.label_map, ...
        seg_2_body_map_temp);
      
      seg_map_auto_2_gt_1 = assign_auto_seg_label_to_groundtruth(...
        seg_auto_1.label_map, body_map_gt_1, map_param);
      seg_map_auto_2_gt_2 = assign_auto_seg_label_to_groundtruth(...
        seg_auto_2.label_map, body_map_gt_2, map_param);
      
      segment_offset = max(linkage_graph_auto.links_3D_p(:,1)) + 1;
      linkage_graph_auto.links_3D_p(:,2) = ...
        linkage_graph_auto.links_3D_p(:,2) + segment_offset;
      seg_map_auto_2_gt_2(:,1) = seg_map_auto_2_gt_2(:,1) + segment_offset;
      
      linkage_graph_auto.links_3D_p = linkage_graph_auto.links_3D_p(...
        apply_mapping(linkage_graph_auto.links_3D_p(:,1), seg_map_auto_2_gt_1)~=0 & ...
        apply_mapping(linkage_graph_auto.links_3D_p(:,2), seg_map_auto_2_gt_2)~=0,:);
      lg = linkage_graph_auto.links_3D_p;
      if(~isempty(lg))
        lg(:,3) = apply_mapping(lg(:,1), seg_map_auto_2_gt_1)==...
          apply_mapping(lg(:,2), seg_map_auto_2_gt_2);
        sec_seg_2_body_map = get_sec_seg_2_body_map_from_links_3D_with_dummy(...
          {lg}, 0.5);
        seg_2_body_map_gt = sec_seg_2_body_map(:, 2:3);
      else
        if(evaluate_config.is_verbose)
          fprintf('seg_2_body_map_gt is empty, skipping\n');
        end
        continue;
      end
      
      if(evaluate_config.is_verbose)
        fprintf('Computing automatic body mapping for linkage thresholds\n');
      end
      rand_score_false_split = zeros(length(evaluate_config.linkage_thresholds),1);
      rand_score_false_split_bound = zeros(length(evaluate_config.linkage_thresholds),1);
      rand_score_false_merge = zeros(length(evaluate_config.linkage_thresholds),1);
      rand_score_false_merge_bound = zeros(length(evaluate_config.linkage_thresholds),1);
      edit_score_n_split = zeros(length(evaluate_config.linkage_thresholds),1);
      edit_score_n_merge = zeros(length(evaluate_config.linkage_thresholds),1);
      n_auto_segment = zeros(length(evaluate_config.linkage_thresholds),1);
      for t = 1:length(evaluate_config.linkage_thresholds)
        sec_seg_2_body_map = get_sec_seg_2_body_map_from_links_3D_with_dummy(...
          {linkage_graph_auto.links_3D_p}, evaluate_config.linkage_thresholds(t));
        seg_2_body_map_auto = sec_seg_2_body_map(:,2:3);
        
        % find automatic segments that are not in the linkage graph and each to a
        % distinct body
        if(~isempty(seg_2_body_map_auto))
          absent_segment_id = unique(seg_2_body_map_gt(...
            ~ismember(seg_2_body_map_gt(:,1), seg_2_body_map_auto(:,1)), 1));
          absent_segment_body_offset = max(seg_2_body_map_auto(:,2))+1;
        else
          absent_segment_id = unique(seg_2_body_map_gt(:,1));
          absent_segment_body_offset = 1;
        end
        seg_2_body_map_auto = [seg_2_body_map_auto; absent_segment_id, ...
          (absent_segment_body_offset:absent_segment_body_offset+length(absent_segment_id)-1)']; %#ok<AGROW>
        
        if(evaluate_config.is_verbose)
          fprintf('Computing rand score\n');
        end
        rand_score = compute_rand_score(uint32(seg_2_body_map_gt), ...
          uint32(seg_2_body_map_auto));
        rand_score_false_split(t) = rand_score(1);
        rand_score_false_split_bound(t) = rand_score(2);
        rand_score_false_merge(t) = rand_score(3);
        rand_score_false_merge_bound(t) = rand_score(4);
        if(evaluate_config.is_verbose)
          fprintf(['linkage_threshold: %g, rand_score_false_split: %d, ' ...
          'rand_score_false_split_bound: %d, rand_score_false_merge: %d, ', ...
          'rand_score_false_merge_bound: %d\n'], evaluate_config.linkage_thresholds(t), ...
          rand_score);
        end
        
        if(evaluate_config.is_verbose)
          fprintf('Computing edit score: #splits and #mergers needed for correction\n');
        end
        [edit_score_1, edit_score_2, split_segment_id, merge_segment_id] = ...
          compute_split_merge_error(seg_2_body_map_gt, seg_2_body_map_auto);
        edit_score_n_split(t) = edit_score_1;
        edit_score_n_merge(t) = edit_score_2;
        n_auto_segment(t) = size(seg_2_body_map_auto,1);
        if(evaluate_config.is_verbose)
          fprintf(['linkage_threshold: %g, edit_score_n_split: %d, ' ...
          'edit_score_n_merge: %d, n_auto_segment: %d\n'], ...
          evaluate_config.linkage_thresholds(t), edit_score_1, edit_score_2, ...
          n_auto_segment(t));
        end
        
        if(evaluate_config.is_verbose_figures)
          im1 = imread([get_stack_dir(config), image_prefixes_1{tile_1}, '.tif']);
          figure(1);
          d = repmat(im1, [1 1 3]);
          d(:,:,1) = uint8(255*ismember(seg_auto_1.label_map, split_segment_id));
          imshow(imresize(d, 1/4));
          title('tile 1 to-be-split segments');
          figure(2);
          d = repmat(im1, [1 1 3]);
          d(:,:,1) = uint8(255*ismember(seg_auto_1.label_map, merge_segment_id));
          imshow(imresize(d, 1/4));
          title('tile 1 to-be-merged segments');
          
          s = seg_auto_2.label_map;
          s(s>0) = s(s>0) + segment_offset;
          im2 = imread([get_stack_dir(config), image_prefixes_2{tile_2}, '.tif']);
          figure(3);
          d = repmat(im2, [1 1 3]);
          d(:,:,1) = uint8(255*ismember(s, split_segment_id));
          imshow(imresize(d, 1/4));
          title('tile 2 to-be-split segments');
          figure(4);
          d = repmat(im2, [1 1 3]);
          d(:,:,1) = uint8(255*ismember(s, merge_segment_id));
          imshow(imresize(d, 1/4));
          title('tile 2 to-be-merged segments');
          
          plot_segment_boundaries(im1, remove_merged_boundaries_2D(uint32(apply_mapping(seg_auto_1.label_map, seg_2_body_map_auto))), 4);
          title('tile 1 auto body map');
          plot_segment_boundaries(im1, remove_merged_boundaries_2D(uint32(body_map_gt_1)), 4);
          title('tile 1 ground-truth body map');
          plot_segment_boundaries(im2, remove_merged_boundaries_2D(uint32(apply_mapping(s, seg_2_body_map_auto))), 4);
          title('tile 2 auto body map');
          plot_segment_boundaries(im2, remove_merged_boundaries_2D(uint32(body_map_gt_2)), 4);
          title('tile 2 ground-truth body map');
        end
      end
      
      link_eval_p.linkage_thresholds = evaluate_config.linkage_thresholds;
      link_eval_p.rand_score_false_split = rand_score_false_split;
      link_eval_p.rand_score_false_split_bound = rand_score_false_split_bound;
      link_eval_p.rand_score_false_merge = rand_score_false_merge;
      link_eval_p.rand_score_false_merge_bound = rand_score_false_merge_bound;
      link_eval_p.edit_score_n_split = edit_score_n_split;
      link_eval_p.edit_score_n_merge = edit_score_n_merge;
      link_eval_p.n_auto_segment = n_auto_segment;
      
      save2(file_name_linkage_eval, 'link_eval_p');
    end
  end
end

if(evaluate_config.is_verbose)
  fprintf('STOP: linkage_evaluate_segment_rand_error_multi_tile\n');
end
return
end
