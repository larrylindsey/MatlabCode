function gen_link_gph_overlap_hist_boundary_intr_hist_sdp_prealign_sift(config)
% gen_link_gph_overlap_hist_boundary_intr_hist_sdp_prealign_sift(config)
% generate linkage graph using optimization on linkage and segmentation
% confidences computed using boost classifier.
% Section alignment :prealigned, single-tile imod, or multi-tile SIFT
% based - see align_stack_SIFT_affine_tile_pair_inter_plane.m.
%
% - Takes as input a stack of images and 2D superpixel and segment maps.
% - Constructs a linkage graph, links_3D{}, with segments as nodes and 3D link
% confidences as the edge weights. This graph's link-weights are
% thresholded in proofread_linkage_superpixel.m to get the final 3D
% segmentation.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  08012009  init code
%

fprintf('START: gen_link_gph_overlap_hist_boundary_intr_hist_sdp_prealign_sift\n');

% Set required paths
p = {
     '~/research/lib/cvx:', ...
     '~/research/lib/cvx/structures:', ...
     '~/research/lib/cvx/lib:', ...
     '~/research/lib/cvx/functions:', ...
     '~/research/lib/cvx/commands:', ...
     '~/research/lib/cvx/builtins:', ...
     '~/research/lib/SDPT3-4.0-beta', ...
     '~/research/lib/SDPT3-4.0-beta/Examples', ...
     '~/research/lib/SDPT3-4.0-beta/HSDSolver', ...
     '~/research/lib/SDPT3-4.0-beta/HSDSolver/etc', ...
     '~/research/lib/SDPT3-4.0-beta/Linsysolver', ...
     '~/research/lib/SDPT3-4.0-beta/Linsysolver/MA47', ...
     '~/research/lib/SDPT3-4.0-beta/Linsysolver/spchol', ...
     '~/research/lib/SDPT3-4.0-beta/Solver', ...
     '~/research/lib/SDPT3-4.0-beta/Solver/Mexfun', ...
     '~/research/lib/SDPT3-4.0-beta/dimacs', ...
     '~/research/lib/SDPT3-4.0-beta/sdplib', ...
     '~/research/lib/SDPT3-4.0-beta/testdir', ...
     '~/research/lib/yalmip', ...
     '~/research/lib/yalmip/demos', ...
     '~/research/lib/yalmip/extras', ...
     '~/research/lib/yalmip/modules', ...
     '~/research/lib/yalmip/modules/global', ...
     '~/research/lib/yalmip/modules/moment', ...
     '~/research/lib/yalmip/modules/parametric', ...
     '~/research/lib/yalmip/modules/robust', ...
     '~/research/lib/yalmip/modules/sos', ...
     '~/research/lib/yalmip/operators', ...
     '~/research/lib/yalmip/solvers', ...
     '~/research/lib/yalmip/usertest', ...
};
for i = 1:length(p)
  addpath(p{i});
end

stack_config = config.stack;
seg_config = config.segmentation_2D(end);
linkage_config = config.linkage;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
apply_config = linkage_config.apply;

version_name = 'overlap_hist_boundary_interior_hist';
version_id = strcmp(feature_config.type, version_name);
if(version_id~=1)
  error('Feature type does not match with called function. Exiting');
end;
if(strcmp(model_config.type, 'boost_sdp')==0)
  error('Model type does not match with called function. Exiting');
end;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(linkage_config.dir, 'dir')~=7)
  mkdir2(linkage_config.dir);
end;
cd(prev_dir);

image_dir = get_stack_dir(config);
sift_dir = get_sift_dir(config);

junk = tree_node_w(1); %#ok<NASGU>
classifier_model = load2([get_reconstruction_dir(config), linkage_config.dir, 'link_model', '.', ...
  model_config.type, '.', feature_config.type, apply_config.model_suffix, '.mat'], ...
  'classifier', 'feature_config', 'annotation_file', 'model_config');

case_id = stack_config.case_ids(1);
fprintf('Constructing linkage graph (links_3D) ..\n%d ', case_id);

% load alignment xf data and save for future
if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned)
  align_config = stack_config.align;
  align_xf_file_name = [image_dir, align_config.xf_file_name];
  xf = read_xf(align_xf_file_name);
  save2([image_dir, align_config.xf_save_name], 'xf');
end;

[images_1, image_prefixes_1] = get_image_from_stack(config, case_id);
if(~isempty(stack_config.roi))
  for tile = 1:length(images_1)
    images_1{tile} = images_1{tile}(stack_config.roi.ymin:stack_config.roi.ymax, ...
      stack_config.roi.xmin:stack_config.roi.xmax);
  end
end;
segment_2D_label_map_1 = {};
for tile = 1:length(images_1)
  images_1{tile} = im2double(images_1{tile});
  [seg_method, seg_suffix] = get_segmentation_suffix(config, image_prefixes_1{tile});
  seg_dir = [get_reconstruction_dir(config), seg_config.dir, seg_method, '/'];
  segment_2D_label_map_1{tile} = load2([seg_dir, image_prefixes_1{tile}, ...
    seg_suffix,'.mat']); %#ok<AGROW>
  if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
    segment_2D_label_map_1{tile}.label_map = imresize(segment_2D_label_map_1{tile}.label_map, ...
      size(images_1{tile}), 'nearest'); %#ok<AGROW>
  end
end

for i = 2:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  case_id_prev = stack_config.case_ids(i-1);
  fprintf('%d:', case_id);
  if(mod(i,20)==0)
    fprintf('\n');
  end;

  [images_2, image_prefixes_2] = get_image_from_stack(config, case_id);
  if(~isempty(stack_config.roi))
    for tile = 1:length(images_2)
      images_2{tile} = images_2{tile}(stack_config.roi.ymin:stack_config.roi.ymax, ...
        stack_config.roi.xmin:stack_config.roi.xmax);
    end
  end;
  segment_2D_label_map_2 = {};
  for tile = 1:length(images_2)
    images_2{tile} = im2double(images_2{tile});
    [seg_method, seg_suffix] = get_segmentation_suffix(config, image_prefixes_2{tile});
    seg_dir = [get_reconstruction_dir(config), seg_config.dir, seg_method, '/'];
    segment_2D_label_map_2{tile} = load2([seg_dir, image_prefixes_2{tile}, ...
      seg_suffix,'.mat']); %#ok<AGROW>
    if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
      segment_2D_label_map_2{tile}.label_map = ...
        imresize(segment_2D_label_map_2{tile}.label_map, ...
        size(images_2{tile}), 'nearest'); %#ok<AGROW>
    end
  end

  links_3D_p = {};

  % compute a linkage graph for each pair of overlapping tiles
  for tile_1 = 1:length(images_1)
    for tile_2 = 1:length(images_2)
      fprintf('%d-%d ', tile_1, tile_2);
      %%%%%%%%%%%%%%%%%%%%
      % For aligning the two images and segmentation maps
      %%%%%%%%%%%%%%%%%%%%
      if(isfield(stack_config, 'image_structure') && ...
          ~isempty(stack_config.image_structure) ) % multi tile system
        file_name_suffix = get_file_name_from_tuple(sift_dir, ...
          image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'st.');
        file_name = [file_name_suffix, '.mat'];
        transforms_tp = [];
        try
          load2(file_name,'transforms_tp');
        catch %#ok<CTCH>
          continue;
        end
        if(isempty(transforms_tp))
          continue;
        end
        xdata = {};
        ydata = {};
        tform_1 = maketform('affine', reshape(transforms_tp(1,:), [2 3])');
        [image_1_tt, xdata{1}, ydata{1}] = ...
          imtransform(images_1{tile_1}, tform_1); %#ok<AGROW>
        minx = min(xdata{1}); miny = min(ydata{1});
        maxx = max(xdata{1}); maxy = max(ydata{1});

        tform_2 = maketform('affine', ...
          reshape(transforms_tp(2,:), [2 3])');
        [image_2_tt, xdata{2}, ydata{2}] = ...
          imtransform(images_2{tile_2}, tform_2); %#ok<AGROW>
        minx = min([minx, xdata{2}]);
        miny = min([miny, ydata{2}]);
        maxx = max([maxx, xdata{2}]);
        maxy = max([maxy, ydata{2}]);
        for t = 1:2
          xdata{t} = round(xdata{t}-minx+1); %#ok<AGROW>
          ydata{t} = round(ydata{t}-miny+1); %#ok<AGROW>
        end
        maxx = round(maxx - minx + 1);
        maxy = round(maxy - miny + 1);

        image_1_t = zeros(maxy, maxx);
        image_1_t(ydata{1}(1):ydata{1}(2), xdata{1}(1):xdata{1}(2)) = ...
          image_1_tt;
        seg_1_tt = imtransform(segment_2D_label_map_1{tile_1}.label_map, tform_1, ...
          'nearest', 'FillValues', -1);
        segment_2D_label_map_1{tile_1}.label_map_t = zeros(maxy, maxx)-1; %#ok<AGROW>
        segment_2D_label_map_1{tile_1}.label_map_t...
          (ydata{1}(1):ydata{1}(2), xdata{1}(1):xdata{1}(2)) = seg_1_tt; %#ok<AGROW>

        image_2_t = zeros(maxy, maxx);
        image_2_t(ydata{2}(1):ydata{2}(2), xdata{2}(1):xdata{2}(2)) = ...
          image_2_tt;
        seg_2_tt = imtransform(segment_2D_label_map_2{tile_2}.label_map, tform_2, ...
          'nearest', 'FillValues', -1);
        segment_2D_label_map_2{tile_2}.label_map_t = zeros(maxy, maxx)-1; %#ok<AGROW>
        segment_2D_label_map_2{tile_2}.label_map_t...
          (ydata{2}(1):ydata{2}(2), xdata{2}(1):xdata{2}(2)) = seg_2_tt; %#ok<AGROW>
      else
        % apply alignment transformation if stack is not pre-aligned
        if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned)
          align_tform = xf_2_tform(xf(i-1,:), size(images_1{tile_1}, 1), ...
            size(images_1{tile_1}, 2));
          margin.left = 1-align_config.margin;
          margin.right = size(images_1{tile_1},2)+align_config.margin;
          margin.top = 1-align_config.margin;
          margin.bottom = size(images_1{tile_1},1)+align_config.margin;

          image_1_t = apply_tform(images_1{tile_1}, align_tform, margin, 'bilinear');
          segment_2D_label_map_1{tile_1}.label_map_t = ...
            apply_tform(segment_2D_label_map_1{tile_1}.label_map, ...
            align_tform, margin, 'nearest'); %#ok<AGROW>

          align_tform = xf_2_tform(xf(i,:), size(images_2{tile_2}, 1), ...
            size(images_2{tile_2}, 2));
          image_2_t = apply_tform(images_2{tile_2}, align_tform, margin, 'bilinear');
          segment_2D_label_map_2{tile_2}.label_map_t = ...
            apply_tform(segment_2D_label_map_2{tile_2}.label_map, ...
            align_tform, margin, 'nearest'); %#ok<AGROW>

          if(isfield(apply_config, 'is_verbose') && apply_config.is_verbose)
            figure(898); imshow(image_1_t);
            figure(899); imshow(image_2_t);
            pause;
          end
        else
          image_1_t = images_1{tile_1};
          image_2_t = images_2{tile_2};
          segment_2D_label_map_1{tile_1}.label_map_t = segment_2D_label_map_1{tile_1}.label_map; %#ok<AGROW>
          segment_2D_label_map_2{tile_2}.label_map_t = segment_2D_label_map_2{tile_2}.label_map; %#ok<AGROW>
        end
      end
      

      %%%
      % build linkage graph for tile pair
      %%%
      links_3D_p{tile_1, tile_2} = []; %#ok<AGROW>
      
      fprintf('relabel the segmentation maps to eliminate non-existent labels and\n');
      fprintf('make labels unique.\n');
      [seg_0, relabelling_0] = relabel_to_remove_nonexistent_labels(...
        segment_2D_label_map_1{tile_1}.label_map_t);
      [seg_1, relabelling_1] = relabel_to_remove_nonexistent_labels(...
        segment_2D_label_map_2{tile_2}.label_map_t);
      relabelling_1(2:end) = relabelling_1(2:end) + max(relabelling_0);
      seg_1(seg_1>0) = seg_1(seg_1>0) + max(relabelling_0);
      
      fprintf('Compute across section linkage confidences\n');
      interior_hist_t = collect_segment_stats_interior_hist(...
        uint32(seg_0), uint8(255*image_1_t), ...
        uint8(feature_config.interior_hist.bin_size));
      interior_hist = [];
      interior_hist(interior_hist_t(:,1)+1, :) = interior_hist_t(:, 2:end); %#ok<AGROW>
      interior_hist_0 = interval_sum(interior_hist);

      interior_hist_t = collect_segment_stats_interior_hist(...
        uint32(seg_1), uint8(255*image_2_t), ...
        uint8(feature_config.interior_hist.bin_size));
      interior_hist = [];
      interior_hist(interior_hist_t(:,1)+1, :) = interior_hist_t(:, 2:end); %#ok<AGROW>
      interior_hist_1 = interval_sum(interior_hist);

      overlap_hist_t = collect_segment_overlap_stats_interior_hist(...
        uint32(seg_0), uint32(seg_1), uint8(255*(image_1_t + image_2_t)/2), ...
        uint8(feature_config.interior_hist.bin_size));
      overlap_hist_t = ...
        overlap_hist_t(overlap_hist_t(:,1)>0 & overlap_hist_t(:,2)>0, :);
      overlap_hist = interval_sum(overlap_hist_t(:, 3:end));

      features_linkage = [overlap_hist, interior_hist_0(1+overlap_hist_t(:,1), :), ...
        interior_hist_1(1+overlap_hist_t(:,2), :)];
      if(size(features_linkage,2)~=classifier_model.feature_config.overlap_hist.n_dimension)
        error('Linkage feature dimensions don''t match during training and testing.');
      end
      weak_learner = tree_node_w(classifier_model.model_config.overlap_hist.tree_depth); %#ok<NASGU>
      confidence_linkage = ...
        Classify(classifier_model.classifier.overlap_hist.GLearners, ...
        classifier_model.classifier.overlap_hist.GWeights, ...
        features_linkage');
      
      fprintf('Compute within section merger confidences\n');
      boundary_hist_0 = collect_segment_pair_stats_boundary_hist(...
        uint32(seg_0), uint8(255*image_1_t), ...
        uint8(classifier_model.feature_config.boundary_hist.bin_size));
      boundary_hist_0 = boundary_hist_0(min(boundary_hist_0(:,1:2), [], 2)>0, :);
      boundary_hist_0 = boundary_hist_0(...
        boundary_hist_0(:,1)~=boundary_hist_0(:,2), :);
      features_boundary = interval_sum(boundary_hist_0(:, 3:end));

      interior_hist_t = collect_segment_stats_interior_hist(...
        uint32(seg_0), uint8(255*image_1_t), ...
        uint8(feature_config.interior_hist.bin_size));
      interior_hist = [];
      interior_hist(interior_hist_t(:,1)+1, :) = interior_hist_t(:, 2:end); %#ok<AGROW>
      interior_hist = interval_sum(interior_hist);
      features_interior = [interior_hist(boundary_hist_0(:,1)+1, :), ...
        interior_hist(boundary_hist_0(:,2)+1, :)];
      
      features = [features_boundary, features_interior];
      if(size(features,2)~=classifier_model.feature_config.boundary_hist.n_dimension)
        error('Segment merge feature dimensions don''t match during training and testing.');
      end
      weak_learner = tree_node_w(classifier_model.model_config.boundary_hist.tree_depth); %#ok<NASGU>
      confidence_boundary_0 = ...
        Classify(classifier_model.classifier.boundary_hist.GLearners, ...
        classifier_model.classifier.boundary_hist.GWeights, ...
        features');

      boundary_hist_1 = collect_segment_pair_stats_boundary_hist(...
        uint32(seg_1), uint8(255*image_2_t), ...
        uint8(classifier_model.feature_config.boundary_hist.bin_size));
      boundary_hist_1 = boundary_hist_1(min(boundary_hist_1(:,1:2), [], 2)>0, :);
      boundary_hist_1 = boundary_hist_1(...
        boundary_hist_1(:,1)~=boundary_hist_1(:,2), :);
      features_boundary = interval_sum(boundary_hist_1(:, 3:end));

      interior_hist_t = collect_segment_stats_interior_hist(...
        uint32(seg_1), uint8(255*image_2_t), ...
        uint8(feature_config.interior_hist.bin_size));
      interior_hist = [];
      interior_hist(interior_hist_t(:,1)+1, :) = interior_hist_t(:, 2:end); %#ok<AGROW>
      interior_hist = interval_sum(interior_hist);
      features_interior = [interior_hist(boundary_hist_1(:,1)+1, :), ...
        interior_hist(boundary_hist_1(:,2)+1, :)];
      
      features = [features_boundary, features_interior];
      if(size(features,2)~=classifier_model.feature_config.boundary_hist.n_dimension)
        error('Segment merge feature dimensions don''t match during training and testing.');
      end
      confidence_boundary_1 = ...
        Classify(classifier_model.classifier.boundary_hist.GLearners, ...
        classifier_model.classifier.boundary_hist.GWeights, ...
        features');

      relabelling_backref_0 = [];
      relabelling_backref_0(relabelling_0(relabelling_0>0)) = ...
        find(relabelling_0)-1; %#ok<AGROW>
      relabelling_backref_1 = [];
      relabelling_backref_1(relabelling_1(relabelling_1>0)) = ...
        find(relabelling_1)-1; %#ok<AGROW>

      min_node_id_0 = 1;
      max_node_id_0 = max(relabelling_0);
      min_node_id_1 = max_node_id_0 + 1;
      max_node_id_1 = max(relabelling_1);
      fprintf('min_node_id_0: %d, max_node_id_0: %d\n', ...
        min_node_id_0, max_node_id_0);
      fprintf('min_node_id_1: %d, max_node_id_1: %d\n', ...
        min_node_id_1, max_node_id_1);
      
      n_node = max_node_id_1;
      fprintf('n_node: %d\n', n_node);
      
      fprintf('build F matrix - symmetric differences\n');
      F = eye(n_node, n_node);
      
      F(sub2ind(size(F), boundary_hist_0(:,1), boundary_hist_0(:,2))) = ...
        -confidence_boundary_0';
      F(sub2ind(size(F), boundary_hist_0(:,2), boundary_hist_0(:,1))) = ...
        -confidence_boundary_0';
      
      F(sub2ind(size(F), boundary_hist_1(:,1), boundary_hist_1(:,2))) = ...
        -confidence_boundary_1';
      F(sub2ind(size(F), boundary_hist_1(:,2), boundary_hist_1(:,1))) = ...
        -confidence_boundary_1';
      
      F(sub2ind(size(F), overlap_hist_t(:,1), ...
        overlap_hist_t(:,2))) = -confidence_linkage';
      F(sub2ind(size(F), overlap_hist_t(:,2), ...
        overlap_hist_t(:,1))) = -confidence_linkage';

      lambda =  1;
      fprintf('lambda: %d\n', lambda);
      Lambda = zeros(size(F))+1;
      Lambda(min_node_id_0:max_node_id_0, min_node_id_0:max_node_id_0) = lambda;
      Lambda(min_node_id_1:max_node_id_1, min_node_id_1:max_node_id_1) = lambda;
      
      G = -(F*100 + Lambda);
      G(G==0) = -1;

%       P = co_cluster_SDP(-G);
      P = co_cluster_SDP_cvx(-G);
      P = post_process_co_cluster_SDP(P);
      
      [node_i, node_j] = find(P);
      valid_flag = node_i>=min_node_id_0 & node_i<=max_node_id_0 & ...
        node_j>=min_node_id_1 & node_j<=max_node_id_1;
      node_0 = node_i(valid_flag);
      node_1 = node_j(valid_flag);
      
      links_3D_p{tile_1, tile_2} = [relabelling_backref_0(node_0(:))', ...
        relabelling_backref_1(node_1(:))', ...
        P(sub2ind(size(P), node_0(:), node_1(:)))]; %#ok<AGROW>
    end
  end
  
  fprintf('Saving links_3D_p for this pair of sections.\n');
  save2([get_reconstruction_dir(config), linkage_config.dir, 'links_3D_p', '.', ...
    num2str(case_id_prev), '.', model_config.type, '.', feature_config.type, ...
    apply_config.model_suffix, config.segmentation_choose.choice.seg_suffix, ...
    '.mat'], 'links_3D_p');

  images_1 = images_2;
  segment_2D_label_map_1 = segment_2D_label_map_2;
  fprintf('\n');
end;

fprintf('done\n');
fprintf('STOP: gen_link_gph_overlap_hist_boundary_intr_hist_sdp_prealign_sift\n');

return;
end
