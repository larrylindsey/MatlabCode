function linkage_3D_train_normalize_overlap_area_boost(config)
% linkage_3D_train_intensity_pair_v2_boost(config)
% Trains a boosted classifier to perform 3D linkage.
%
% - Features are collected for spatially overlapping 2D segments in manually
% annotated data. If two segments have the same label then their feature is
% put in the postive training set, else in the negative training set.
% - Currently the features are area of overlap
%
% Pratim Ghosh,
% Summer intern, Janelia Farm Research Campus, HHMI.
% Univ. of California, Santa Barbara.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  06202009 init code
%

fprintf('START: linkage_3D_train_normalize_overlap_area_boost\n');
linkage_config = config.linkage;
train_config = linkage_config.train;
feature_config = linkage_config.feature;
model_config = linkage_config.model;

if(strcmp(feature_config.type, 'normalize_overlap_area')==0)
  error('Feature type does not match with called function. Exiting');
end;
if(strcmp(model_config.type, 'boost')==0)
  error('Classifier type does not match with called function. Exiting');
end;

% Check to make sure that linkage directory exists.
prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(linkage_config.dir, 'dir')~=7)
  mkdir2(linkage_config.dir);
end;
cd(prev_dir);

fprintf('Loading manual annotation .. ');
manual_annotation = load2(train_config.manual_annotation_file);
fprintf('done\n');

fprintf('Collecting linkage statistics .. ');
train_labels = [];
features = [];
for z = 1:length(manual_annotation.al)-1
  %%%%%%%%%%%%%%%%%%%%%%%% CHANGE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Task: Compute features between pairs of overlapping segments in section
  % z and z+1.
  % To generate:
  % features{z}: feature vectors, one for each pair of overlapping
  % segments. A NxD matrix where N is the number of overlapping pairs of
  % segments, and D is the dimensionality of the feature.
  %
  % Given:
  % (a) manual_annotation.al{z}: image for section z
  % (b) manual_annotation.wcat{z}: body label map for section z
  % (c) manual_annotation.al{z+1}: image for section z+1
  % (d) manual_annotation.wcat{z+1}: body label map for section z+1
  % (e) unique_label_pairs: pairs of labels of overlapping segments. Should
  % correspond to features{z}. A Nx2 matrix [labels_1, labels_2].
  % features{z}(i,:) should the feature vector for unique_label_pairs(i,:)
  
  % features
  label_pairs = [reshape(manual_annotation.wcat{z}, [numel(manual_annotation.al{z}), 1]), ...
    reshape(manual_annotation.wcat{z+1}, [numel(manual_annotation.al{z}), 1])];
  label_pairs = label_pairs(label_pairs(:,1)>0 & label_pairs(:,2)>0, :);
  
  sorted_label_pairs =sortrows(label_pairs); %% added
  [unique_label_pairs I]=unique(sorted_label_pairs,'rows');
  [junk I_f]=unique(sorted_label_pairs,'rows','first');
  overlap_area=I-I_f+1;
  
  
  label1=unique_label_pairs(:,1);
  label2=unique_label_pairs(:,2);
  
  unique_label1=unique(label1);
  unique_label2=unique(label2);
  
  junk1=hist(sorted_label_pairs(:,1),1:max(unique_label1));
  junk2=hist(sorted_label_pairs(:,2),1:max(unique_label2));
  
  individual_area1=(junk1(label1))';
  individual_area2=(junk2(label2))';
  
   
  features{z} = overlap_area ./ sqrt(individual_area1.*individual_area2); %#ok<AGROW>
  %%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE DONE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % +1 if segments have same labels and -1 otherwise.
  train_labels{z} = ...
    2*((unique_label_pairs(:,1)==unique_label_pairs(:,2)) - 0.5); %#ok<AGROW>
end; 

fprintf('done\n');

fprintf('Training classifier .. ');
features_trn = [];
labels_trn = [];
for i = 1:length(features)
  features_trn = [features_trn; features{i}]; %#ok<AGROW>
  labels_trn = [labels_trn; train_labels{i}]; %#ok<AGROW>
end;

MaxIter = model_config.n_iteration;
tree_depth = model_config.tree_depth;
weak_learner = tree_node_w(tree_depth);
[GLearners, GWeights] = ...
  GentleAdaBoost(weak_learner, features_trn', labels_trn', MaxIter); %#ok<NASGU>

annotation_file = train_config.manual_annotation_file; %#ok<NASGU>
save2([get_reconstruction_dir(config), linkage_config.dir, 'link_model', '.', ...
  model_config.type, '.', feature_config.type, train_config.save_suffix, '.mat'], ...
  'GLearners', 'GWeights', 'MaxIter', 'tree_depth', 'annotation_file');
fprintf('done\n');

fprintf('STOP: linkage_3D_train_normalize_overlap_area_boost\n');
return
end
