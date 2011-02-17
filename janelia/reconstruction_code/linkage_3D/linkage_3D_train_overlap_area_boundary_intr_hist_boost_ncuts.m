function linkage_3D_train_overlap_area_boundary_intr_hist_boost_ncuts(config)
% linkage_3D_train_overlap_area_boundary_hist_boost_lp_p3(config)
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

fprintf('START: linkage_3D_train_overlap_area_boundary_intr_hist_boost_ncuts\n');

linkage_config = config.linkage;
train_config = linkage_config.train;
feature_config = linkage_config.feature;
model_config = linkage_config.model;

version_names{1} = 'overlap_area_boundary_interior_hist';
version_id = strcmp(feature_config.type, version_names{1});
if(version_id~=1)
  error('Feature type does not match with called function. Exiting');
end;
if(strcmp(model_config.type, 'boost_ncuts')==0)
  error('Classifier type does not match with called function. Exiting');
end;

linkage_dir = get_linkage_dir(config);

fprintf('Loading manual annotation .. ');
manual_annotation = load2(train_config.manual_annotation_file);
fprintf('done\n');

% for each pair of adjacent slices in the manual annotations
%   for each pair of overlapping segments in the adjacent slices
%     - collect feature: <overlap area, area1, area2>
%     - if the two segments belong to the same process then put feature in
%     postive set.
%     - else put the feature in negative set.
fprintf('Collecting linkage statistics ..\n');
labels = [];
features = [];
for z = 1:length(manual_annotation.al)-1
  fprintf('z: %d\n', z);
  label_pairs = [reshape(manual_annotation.wcat{z}, [numel(manual_annotation.al{z}), 1]), ...
    reshape(manual_annotation.wcat{z+1}, [numel(manual_annotation.al{z}), 1])];
  label_pairs = label_pairs(label_pairs(:,1)>0 & label_pairs(:,2)>0, :);
  [label_pairs_unique, label_pairs_count] = count_row_occurence(label_pairs);
  
  [s_1, a_1] = ...
    count_row_occurence(manual_annotation.wcat{z}(:));
  areas_1 = [];
  areas_1(s_1+1) = a_1; %#ok<AGROW>
  [s_2, a_2] = ...
    count_row_occurence(manual_annotation.wcat{z+1}(:));
  areas_2 = [];
  areas_2(s_2+1) = a_2; %#ok<AGROW>
  
  features{z} = [areas_1(1+label_pairs_unique(:,1))', ...
    areas_2(1+label_pairs_unique(:,2))', label_pairs_count; ...
    areas_2(1+label_pairs_unique(:,2))', ...
    areas_1(1+label_pairs_unique(:,1))', label_pairs_count]; %#ok<AGROW>
  
  if(isfield(manual_annotation, 'label_mapping_as_to_gt') && ...
      ~isempty(manual_annotation.label_mapping_as_to_gt))
    is_same_label = ...
      manual_annotation.label_mapping_as_to_gt{z}(1+label_pairs_unique(:,1)) ...
      == manual_annotation.label_mapping_as_to_gt{z+1}(1+label_pairs_unique(:,2));
  else
    is_same_label = label_pairs_unique(:,1) == label_pairs_unique(:,2);
  end
  labels{z} = [2*double(is_same_label)-1; 2*double(is_same_label)-1]; %#ok<AGROW>
end
fprintf('done\n');

%
% train a adaboost based classifier
%
fprintf('Training classifier for linkage across sections ..\n');
features_trn = [];
labels_trn = [];
for i = 1:length(features)
  features_trn = [features_trn; features{i}]; %#ok<AGROW>
  labels_trn = [labels_trn; labels{i}]; %#ok<AGROW>
end;
fprintf('done.\n');
relative_error_weights = ones(size(labels_trn));
if(isfield(model_config.overlap_area, 'merge_error_relative_weight'))
  relative_error_weights(labels_trn==-1) = ...
    model_config.overlap_area.merge_error_relative_weight;
end

feature_config.overlap_area.n_dimension = size(features_trn, 2);
fprintf('Collected %d feature vectors of dimensionality %d\n', ...
  size(features_trn,1), size(features_trn,2));

MaxIter = model_config.overlap_area.n_iteration;
tree_depth = model_config.overlap_area.tree_depth; % 25,2
weak_learner = tree_node_w(tree_depth);
[classifier.overlap_area.GLearners, classifier.overlap_area.GWeights] = ...
  GentleAdaBoost(weak_learner, features_trn', labels_trn', MaxIter, ...
  relative_error_weights');

% Train a classifier to compute confidence values for mergers within a
% section. The features are histograms of boundary values between
% segments.

fprintf('Collecting within-section segmentation statistics ..\n');
labels = [];
features = [];
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
  
  features{z} = [[boundary_feature; boundary_feature], interior_feature]; %#ok<AGROW>

  if(isfield(manual_annotation, 'label_mapping_as_to_gt') && ...
      ~isempty(manual_annotation.label_mapping_as_to_gt))
    is_same_label = ...
      manual_annotation.label_mapping_as_to_gt{z}(1+boundary_hist(:,1)) ...
      == manual_annotation.label_mapping_as_to_gt{z}(1+boundary_hist(:,2));
  else
    is_same_label = boundary_hist(:,1) == boundary_hist(:,2);
  end
  labels{z} = [2*double(is_same_label)-1; 2*double(is_same_label)-1]; %#ok<AGROW>
end
fprintf('done\n');

%
% train a adaboost based classifier
%
features_trn = [];
labels_trn = [];
for i = 1:length(features)
  features_trn = [features_trn; features{i}]; %#ok<AGROW>
  labels_trn = [labels_trn; labels{i}]; %#ok<AGROW>
end;
relative_error_weights = ones(size(labels_trn));
if(isfield(model_config.boundary_hist, 'merge_error_relative_weight'))
  relative_error_weights(labels_trn==-1) = ...
    model_config.boundary_hist.merge_error_relative_weight;
end

feature_config.boundary_hist.n_dimension = size(features_trn, 2);
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

fprintf('STOP: linkage_3D_train_overlap_area_boundary_intr_hist_boost_ncuts\n');
return;
end
