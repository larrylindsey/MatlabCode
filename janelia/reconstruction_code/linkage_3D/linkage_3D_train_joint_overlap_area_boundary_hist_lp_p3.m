function linkage_3D_train_joint_overlap_area_boundary_hist_lp_p3(config)
% linkage_3D_train_symm_diff_mean_boundary_lp_p3(config)
% Trains to perform 3D linkage.
% Method: minimize pixelwise symmetric difference between linked segments.
% Optimization performed using linear programming formulation. See
% "Clustering with qualitative constraints", Charikar et al.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  08012009  init code
%

fprintf('START: linkage_3D_train_symm_diff_mean_boundary_lp_p3\n');

linkage_config = config.linkage;
train_config = linkage_config.train;
feature_config = linkage_config.feature;
model_config = linkage_config.model;

version_names{1} = 'intensity_pair_hist_v2';
version_names{2} = 'intensity_pair_hist_v2b';
version_names{3} = 'intensity_pair_hist_v2c';
version_names{4} = 'intensity_pair_hist_v2d';
version_id = find(strcmp(feature_config.type, version_names));
if(isempty(version_id))
  error('Feature type does not match with called function. Exiting');
end;
if(strcmp(model_config.type, 'boost')==0)
  error('Classifier type does not match with called function. Exiting');
end;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(linkage_config.dir, 'dir')~=7)
  mkdir2(linkage_config.dir);
end;
cd(prev_dir);

fprintf('Loading manual annotation .. ');
manual_annotation = load2(train_config.manual_annotation_file);
fprintf('done\n');

% bins used in the intensity histograms
intensity_bins = feature_config.intensity_bins;

% for each pair of adjacent slices in the manual annotations
%   for each pair of overlapping segments in the adjacent slices
%     - collect pairs of intensities in the segments {<I(x,y,z),I(x,y,z+1)>}
%     and compute a 2D histogram for this set of 2D vectors. This histogram
%     is the feature for 3D linkage.
%     - if the two segments belong to the same process then put feature in
%     postive set.
%     - else put the feature in negative set.
fprintf('Collecting linkage statistics .. ');
labels = [];
label_pairs_all = [];
features = [];
for z = 1:length(manual_annotation.al)-1
  image_1 = histeq(im2double(manual_annotation.al{z}));
  image_2 = histeq(im2double(manual_annotation.al{z+1}));
  intensity_pairs = [reshape(image_1, [numel(manual_annotation.al{z}), 1]), reshape(image_2, [numel(manual_annotation.al{z}), 1])];
  label_pairs = [reshape(manual_annotation.wcat{z}, [numel(manual_annotation.al{z}), 1]), ...
    reshape(manual_annotation.wcat{z+1}, [numel(manual_annotation.al{z}), 1])];
  label_props_1 = regionprops(manual_annotation.wcat{z}, 'Area');
  label_props_2 = regionprops(manual_annotation.wcat{z+1}, 'Area');
  intensity_pairs = intensity_pairs(label_pairs(:,1)>0 & label_pairs(:,2)>0, :);
  label_pairs = label_pairs(label_pairs(:,1)>0 & label_pairs(:,2)>0, :);
  
  [unique_label_pairs, I, J] = unique(label_pairs, 'rows');
  label_pairs_all{z} = unique_label_pairs; %#ok<AGROW>
  
  for p = 1:size(unique_label_pairs, 1)
    unique_label_pairs(p, 3) = sum(J==p);
  end;
  
  if(isfield(manual_annotation, 'label_mapping_as_to_gt') && ...
      ~isempty(manual_annotation.label_mapping_as_to_gt))
    is_same_label = ...
      manual_annotation.label_mapping_as_to_gt{z}(1+unique_label_pairs(:,1)) ...
      == manual_annotation.label_mapping_as_to_gt{z+1}(1+unique_label_pairs(:,2));
  else
    is_same_label = unique_label_pairs(:,1) == unique_label_pairs(:,2);
  end
  labels{z} = [2*double(is_same_label)-1; 2*double(is_same_label)-1]; %#ok<AGROW>

  template_hist = ones(length(intensity_bins));
  template_p = pyramid_hist2(template_hist,2);
  switch(version_id)
    case 1
      intensity_pair_hist_feature = zeros(size(unique_label_pairs,1), length(intensity_bins)^2 + length(template_p));
      intensity_pair_hist_feature_t = zeros(size(unique_label_pairs,1), length(intensity_bins)^2  + length(template_p));
    case 2
      intensity_pair_hist_feature = zeros(size(unique_label_pairs,1), length(intensity_bins)^2 + length(template_p) + 3);
      intensity_pair_hist_feature_t = zeros(size(unique_label_pairs,1), length(intensity_bins)^2  + length(template_p) + 3);
    case 3
      intensity_pair_hist_feature = zeros(size(unique_label_pairs,1), 3);
      intensity_pair_hist_feature_t = zeros(size(unique_label_pairs,1), 3);
    case 4
      intensity_pair_hist_feature = zeros(size(unique_label_pairs,1), 5);
      intensity_pair_hist_feature_t = zeros(size(unique_label_pairs,1), 5);
  end
  for p = 1:size(unique_label_pairs, 1)
    intensity_pair_set = intensity_pairs(J==p, :);
    switch(version_id)
      case 1
        h = hist2(intensity_pair_set, intensity_bins, intensity_bins);
        h_p = pyramid_hist2(h,2);
        h_p_t = pyramid_hist2(h',2);
        intensity_pair_hist_feature(p,:) = [reshape(h, [1, numel(h)]), h_p];
        intensity_pair_hist_feature_t(p,:) = [reshape(h', [1, numel(h)]), h_p_t];
      case 2
        h = hist2(intensity_pair_set, intensity_bins, intensity_bins);
        h_p = pyramid_hist2(h,2);
        h_p_t = pyramid_hist2(h',2);
        intensity_pair_hist_feature(p,:) = [label_props_1(unique_label_pairs(p,1)).Area, ...
          label_props_2(unique_label_pairs(p,2)).Area, size(intensity_pair_set, 1), ...
          reshape(h, [1, numel(h)]), h_p];
        intensity_pair_hist_feature_t(p,:) = [label_props_2(unique_label_pairs(p,2)).Area, ...
          label_props_1(unique_label_pairs(p,1)).Area, size(intensity_pair_set, 1), ...
          reshape(h', [1, numel(h)]), h_p_t];
      case 3
        intensity_pair_hist_feature(p,:) = [label_props_1(unique_label_pairs(p,1)).Area, ...
          label_props_2(unique_label_pairs(p,2)).Area, size(intensity_pair_set, 1)];
        intensity_pair_hist_feature_t(p,:) = [label_props_2(unique_label_pairs(p,2)).Area, ...
          label_props_1(unique_label_pairs(p,1)).Area, size(intensity_pair_set, 1)];
      case 4
        min_segment_area = min(label_props_1(unique_label_pairs(p,1)).Area, ...
          label_props_2(unique_label_pairs(p,2)).Area);
        intensity_pair_hist_feature(p,:) = [label_props_1(unique_label_pairs(p,1)).Area, ...
          label_props_2(unique_label_pairs(p,2)).Area, size(intensity_pair_set, 1), ...
          min_segment_area-size(intensity_pair_set, 1), ...
          size(intensity_pair_set, 1)/min_segment_area];
        intensity_pair_hist_feature_t(p,:) = [label_props_2(unique_label_pairs(p,2)).Area, ...
          label_props_1(unique_label_pairs(p,1)).Area, size(intensity_pair_set, 1), ...
          min_segment_area-size(intensity_pair_set, 1), ...
          size(intensity_pair_set, 1)/min_segment_area];
    end
  end;
  features{z} = [intensity_pair_hist_feature; intensity_pair_hist_feature_t]; %#ok<AGROW>
end;
fprintf('done\n');

%
% train a adaboost based classifier
%
fprintf('Training classifier .. ');
features_trn = [];
labels_trn = [];
for i = 1:length(features)
  features_trn = [features_trn; features{i}]; %#ok<AGROW>
  labels_trn = [labels_trn; labels{i}]; %#ok<AGROW>
end;

MaxIter = model_config.n_iteration;
tree_depth = model_config.tree_depth; % 25,2
weak_learner = tree_node_w(tree_depth);
[GLearners, GWeights] = ...
  GentleAdaBoost(weak_learner, features_trn', labels_trn', MaxIter); %#ok<NASGU>

annotation_file = train_config.manual_annotation_file; %#ok<NASGU>
save2([get_reconstruction_dir(config), linkage_config.dir, 'link_model', '.', ...
  model_config.type, '.', feature_config.type, train_config.save_suffix, '.mat'], ...
  'GLearners', 'GWeights', 'MaxIter', 'tree_depth', 'annotation_file', 'intensity_bins');
fprintf('done\n');

fprintf('STOP: linkage_3D_train_symm_diff_mean_boundary_lp_p3\n');
return;
end
