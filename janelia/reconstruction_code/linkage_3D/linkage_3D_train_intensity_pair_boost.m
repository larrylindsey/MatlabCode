function linkage_3D_train_intensity_pair_boost(config)
% linkage_3D_train_intensity_pair_boost(config)
% Trains a boosted classifier to perform 3D linkage.
%
% - Features are collected for spatially overlapping 2D segments in manually
% annotated data. If two segments have the same label then their feature is
% put in the postive training set, else in the negative training set.
% - Currently the features are intensity histograms
% - A boosted classifier is trained on these features. The classifier must
% be saved to mat file for use during 3D linkage.
%
% Dependencies:
%   1. hist2(.) a mex function for computing 2D histograms. See hist2.cpp
%   for details.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  ~03202008 init code
% v1  04112008  modified for reconstruction pipeline
%

linkage_config = config.linkage;
train_config = linkage_config.train;
feature_config = linkage_config.feature;
model_config = linkage_config.model;

if(strcmp(feature_config.type, 'intensity_pair_hist')==0)
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
manual_annotation = load(train_config.manual_annotation_file);
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

  intensity_pairs = intensity_pairs(label_pairs(:,1)>0 & label_pairs(:,2)>0, :);
  label_pairs = label_pairs(label_pairs(:,1)>0 & label_pairs(:,2)>0, :);
  
  [unique_label_pairs, I, J] = unique(label_pairs, 'rows');
  label_pairs_all{z} = unique_label_pairs;
  
  for p = 1:size(unique_label_pairs, 1)
    unique_label_pairs(p, 3) = sum(J==p);
  end;
  
  is_same_label = unique_label_pairs(:,1) == unique_label_pairs(:,2);
  labels{z} = [2*double(is_same_label)-1; 2*double(is_same_label)-1];

  intensity_pair_hist_feature = zeros(size(unique_label_pairs,1), length(intensity_bins)^2);
  intensity_pair_hist_feature_t = zeros(size(unique_label_pairs,1), length(intensity_bins)^2);
  for p = 1:size(unique_label_pairs, 1)
    intensity_pair_set = intensity_pairs(J==p, :);
    h = hist2(intensity_pair_set, intensity_bins, intensity_bins);
    intensity_pair_hist_feature(p,:) = reshape(h, [1, numel(h)]);
    intensity_pair_hist_feature_t(p,:) = reshape(h', [1, numel(h)]);
  end;
  features{z} = [intensity_pair_hist_feature; intensity_pair_hist_feature_t;];
end;
fprintf('done\n');

%
% train a adaboost based classifier
%
fprintf('Training classifier .. ');
features_trn = [];
labels_trn = [];
for i = 1:length(features)
  features_trn = [features_trn; features{i}];
  labels_trn = [labels_trn; labels{i}];
end;

MaxIter = model_config.n_iteration;
tree_depth = model_config.tree_depth; % 25,2
weak_learner = tree_node_w(tree_depth);
[GLearners, GWeights] = GentleAdaBoost(weak_learner, features_trn', labels_trn', MaxIter);

annotation_file = train_config.manual_annotation_file;
save([get_reconstruction_dir(config), linkage_config.dir, 'link_model', '.', ...
  model_config.type, '.', feature_config.type, train_config.save_suffix, '.mat'], ...
  'GLearners', 'GWeights', 'MaxIter', 'tree_depth', 'annotation_file', 'intensity_bins');
fprintf('done\n');
