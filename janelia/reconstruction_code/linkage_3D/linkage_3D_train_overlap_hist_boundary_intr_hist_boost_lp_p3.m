function linkage_3D_train_overlap_hist_boundary_intr_hist_boost_lp_p3(config)
% linkage_3D_train_overlap_hist_boundary_hist_boost_lp_p3(config)
% Trains to perform 3D linkage.
% Method: Train boost classifiers to link segments in adjacent sections
% based on area of overlap and their areas, and segments within a section
% based on the histogram in boundary values.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  08012009  init code
%

fprintf('START: linkage_3D_train_overlap_hist_boundary_intr_hist_boost_lp_p3\n');
stack_config = config.stack;
linkage_config = config.linkage;
train_config = linkage_config.train;
feature_config = linkage_config.feature;
model_config = linkage_config.model;

version_names{1} = 'overlap_hist_boundary_interior_hist';
version_id = strcmp(feature_config.type, version_names{1});
if(version_id~=1)
  error('Feature type does not match with called function. Exiting');
end;
if(max(strcmp(model_config.type, {'boost_lp_p3', 'boost_sdp'}))==0)
  error('Classifier type does not match with called function. Exiting');
end;

linkage_dir = get_linkage_dir(config);

fprintf('Loading manual annotation .. ');
manual_annotation = load2(train_config.manual_annotation_file);
fprintf('done\n');

if(isfield(train_config, 'is_trained_on_automatic_segmentation') && ...
    train_config.is_trained_on_automatic_segmentation && ...
    (~isfield(manual_annotation, 'label_mapping_as_to_gt') || ...
    isempty(manual_annotation.label_mapping_as_to_gt)))
  fprintf('Loading automatic segmentations for training.\n');
  manual_annotation.wcat_gt = manual_annotation.wcat;
  if(length(stack_config.case_ids)~=length(manual_annotation.wcat_gt))
    error('During training the case_ids should be same as those for ground-truth\n');
  end
  manual_annotation.wcat = {};
  for i = 1:length(stack_config.case_ids)
    case_id = stack_config.case_ids(i);
    fprintf('case_id: %d\n', case_id);
    image_prefixes = get_image_prefixes_subdirs(config, case_id);
    image_prefix = image_prefixes{1};
    [seg_method, seg_suffix] = get_segmentation_suffix(config, image_prefix);
    seg_dir = [get_reconstruction_dir(config), ...
      config.segmentation_2D.dir, seg_method, '/'];
    sg = load2([seg_dir, image_prefix, seg_suffix,'.mat']);
    manual_annotation.wcat{i} = sg.label_map;
  end
  
  fprintf('Assigning automatic segmentation to ground-truth segments ...\n');
  manual_annotation.label_mapping_as_to_gt = {};
  for z = 1:length(stack_config.case_ids)
    fprintf('z: %d\n', z);
    % Assign automatic segmentation map's segments to ground-truth segments.
    % If an A.S. segment overlaps substantially with more than one
    % ground-truth segment then disregard it during training and testing.
    [label_pairs, overlap_area] = count_row_occurence(...
      [manual_annotation.wcat{z}(:), manual_annotation.wcat_gt{z}(:)]);
    overlap_area = overlap_area(label_pairs(:,1)>0 & label_pairs(:,2)>0);
    label_pairs = label_pairs(label_pairs(:,1)>0 & label_pairs(:,2)>0, :);
    
    [junk, a_1] = count_row_occurence(manual_annotation.wcat{z}(:));
    area_1 = [];
    area_1(1+junk) = a_1; %#ok<AGROW>
    [junk, a_2] = count_row_occurence(manual_annotation.wcat_gt{z}(:));
    area_2 = [];
    area_2(1+junk) = a_2; %#ok<AGROW>
    mapping = correspond_labels_many_to_one(label_pairs(:,1), ...
      label_pairs(:,2), overlap_area, area_1', area_2', 200, 0.1);
    manual_annotation.label_mapping_as_to_gt{z} = ...
      zeros(max(manual_annotation.wcat{z}(:))+1, 1);
    manual_annotation.label_mapping_as_to_gt{z}(mapping(:,1) + 1) = ...
      mapping(:,2);
    
    is_mapped_to_0 = find(manual_annotation.label_mapping_as_to_gt{z}==0) - 1;
    manual_annotation.wcat{z}(ismember(manual_annotation.wcat{z}, is_mapped_to_0)) = 0;
  end
  
  manual_annotation.wcat_gt = {};
end

% for each pair of adjacent slices in the manual annotations
%   for each pair of overlapping segments in the adjacent slices
%     - collect feature: histogram of intensity values in the overlapping
%     region and within the two segments.
%     - if the two segments belong to the same process then put feature in
%     postive set.
%     - else put the feature in negative set.
fprintf('Collecting linkage statistics ..\n');
features_trn = [];
labels_trn = [];
for z = 1:length(manual_annotation.al)-1
  fprintf('z: %d\n', z);
  
  interior_hist_t = collect_segment_stats_interior_hist(...
    uint32(manual_annotation.wcat{z}), ...
    uint8(manual_annotation.al{z}), ...
    uint8(feature_config.interior_hist.bin_size));
  interior_hist = [];
  interior_hist(interior_hist_t(:,1)+1, :) = interior_hist_t(:, 2:end); %#ok<AGROW>
  interior_hist_0 = interval_sum(interior_hist);

  interior_hist_t = collect_segment_stats_interior_hist(...
    uint32(manual_annotation.wcat{z+1}), ...
    uint8(manual_annotation.al{z+1}), ...
    uint8(feature_config.interior_hist.bin_size));
  interior_hist = [];
  interior_hist(interior_hist_t(:,1)+1, :) = interior_hist_t(:, 2:end); %#ok<AGROW>
  interior_hist_1 = interval_sum(interior_hist);
  
  overlap_hist_t = collect_segment_overlap_stats_interior_hist(...
    uint32(manual_annotation.wcat{z}), uint32(manual_annotation.wcat{z+1}), ...
    uint8((double(manual_annotation.al{z}) + double(manual_annotation.al{z+1}))/2), ...
    uint8(feature_config.interior_hist.bin_size));
  overlap_hist_t = overlap_hist_t(overlap_hist_t(:,1)>0 & overlap_hist_t(:,2)>0, :);
  overlap_hist = interval_sum(overlap_hist_t(:, 3:end));

  features_trn = [features_trn; ...
    overlap_hist, interior_hist_0(1+overlap_hist_t(:,1),:), ...
    interior_hist_1(1+overlap_hist_t(:,2),:); ...
    overlap_hist, interior_hist_1(1+overlap_hist_t(:,2),:), ...
    interior_hist_0(1+overlap_hist_t(:,1),:)]; %#ok<AGROW>
  
  if(isfield(manual_annotation, 'label_mapping_as_to_gt') && ...
      ~isempty(manual_annotation.label_mapping_as_to_gt))
    is_same_label = ...
      manual_annotation.label_mapping_as_to_gt{z}(1+overlap_hist_t(:,1)) ...
      == manual_annotation.label_mapping_as_to_gt{z+1}(1+overlap_hist_t(:,2));
  else
    is_same_label = overlap_hist_t(:,1) == overlap_hist_t(:,2);
  end
  labels_trn = [labels_trn; ...
    2*double(is_same_label)-1; 2*double(is_same_label)-1]; %#ok<AGROW>
end
fprintf('done\n');

%
% train a adaboost based classifier
%
relative_error_weights = ones(size(labels_trn));
if(isfield(model_config.overlap_hist, 'merge_error_relative_weight'))
  relative_error_weights(labels_trn==-1) = ...
    model_config.overlap_hist.merge_error_relative_weight;
end

feature_config.overlap_hist.n_dimension = size(features_trn,2);
fprintf('Collected %d feature vectors of dimensionality %d\n', ...
  size(features_trn,1), size(features_trn,2));

fprintf('Training classifier for linkage across sections ..\n');
MaxIter = model_config.overlap_hist.n_iteration;
tree_depth = model_config.overlap_hist.tree_depth; % 25,2
weak_learner = tree_node_w(tree_depth);
[classifier.overlap_hist.GLearners, classifier.overlap_hist.GWeights] = ...
  GentleAdaBoost(weak_learner, features_trn', labels_trn', MaxIter, ...
  relative_error_weights');
fprintf('done.\n');

% Train a classifier to compute confidence values for mergers within a
% section. The features are histograms of boundary values between
% segments.
fprintf('Collecting within-section segmentation statistics ..\n');
features_trn = [];
labels_trn = [];
for z = 1:length(manual_annotation.al)
  fprintf('z: %d\n', z);
  boundary_hist = collect_segment_pair_stats_boundary_hist(...
    uint32(manual_annotation.wcat{z}), ...
    uint8(manual_annotation.al{z}), ...
    uint8(feature_config.boundary_hist.bin_size));
  boundary_feature = interval_sum(boundary_hist(:, 3:end));
  
  interior_hist_t = collect_segment_stats_interior_hist(...
    uint32(manual_annotation.wcat{z}), ...
    uint8(manual_annotation.al{z}), ...
    uint8(feature_config.interior_hist.bin_size));
  interior_hist = [];
  interior_hist(interior_hist_t(:,1)+1, :) = interior_hist_t(:, 2:end); %#ok<AGROW>
  interior_hist = interval_sum(interior_hist);
  interior_feature = [interior_hist(boundary_hist(:,1)+1, :), ...
    interior_hist(boundary_hist(:,2)+1, :); ...
    interior_hist(boundary_hist(:,2)+1, :), ...
    interior_hist(boundary_hist(:,1)+1, :)];
  
  features_trn = [features_trn; ...
    [boundary_feature; boundary_feature], interior_feature]; %#ok<AGROW>

  if(isfield(manual_annotation, 'label_mapping_as_to_gt') && ...
      ~isempty(manual_annotation.label_mapping_as_to_gt))
    is_same_label = ...
      manual_annotation.label_mapping_as_to_gt{z}(1+boundary_hist(:,1)) ...
      == manual_annotation.label_mapping_as_to_gt{z}(1+boundary_hist(:,2));
  else
    is_same_label = boundary_hist(:,1) == boundary_hist(:,2);
  end
  labels_trn = [labels_trn; ...
    2*double(is_same_label)-1; 2*double(is_same_label)-1]; %#ok<AGROW>
end
fprintf('done\n');

%
% train a adaboost based classifier
%
relative_error_weights = ones(size(labels_trn));
if(isfield(model_config.boundary_hist, 'merge_error_relative_weight'))
  relative_error_weights(labels_trn==-1) = ...
    model_config.boundary_hist.merge_error_relative_weight;
end
feature_config.boundary_hist.n_dimension = size(features_trn,2);
fprintf('Collected %d feature vectors of dimensionality %d\n', ...
  size(features_trn,1), size(features_trn,2));

fprintf('Training classifier for within-section mergers ..\n');
MaxIter = model_config.boundary_hist.n_iteration;
tree_depth = model_config.boundary_hist.tree_depth; % 25,2
weak_learner = tree_node_w(tree_depth);
[classifier.boundary_hist.GLearners, classifier.boundary_hist.GWeights] = ...
  GentleAdaBoost(weak_learner, features_trn', labels_trn', MaxIter, ...
  relative_error_weights');
fprintf('done.\n');

annotation_file = train_config.manual_annotation_file; %#ok<NASGU>
save2([linkage_dir, 'link_model', '.', ...
  model_config.type, '.', feature_config.type, train_config.save_suffix, '.mat'], ...
  'classifier', 'model_config', 'feature_config', 'annotation_file');
fprintf('done\n');

fprintf('STOP: linkage_3D_train_overlap_hist_boundary_intr_hist_boost_lp_p3\n');
return;
end
