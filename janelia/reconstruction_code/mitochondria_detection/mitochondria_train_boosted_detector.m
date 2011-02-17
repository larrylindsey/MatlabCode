function mitochondria_train_boosted_detector(config)
% mitochondria_train_boosted_detector(config)
% train an boosted classifier for mitochondria detection.
%
% Uses GML Adaboost Library.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  ~Feb.2008  init code
% v1  04112008  modified for reconstruction pipeline
%

MAX_N_VECT = 500000;

is_debugging = config.DEBUG;

mitochondria_config = config.mitochondria;
train_config = mitochondria_config.train;
feature_config = mitochondria_config.feature;
model_config = mitochondria_config.model;

if(strcmp(model_config.type, 'boost')==0)
  error('Classifier type does not match with called function. Exiting');
end;
if(~isfield(model_config, 'version'))
  model_config.version = '';
end

if(~isfield(feature_config, 'version'))
  feature_config.version = '';
end

image_dir = [get_reconstruction_dir(config), mitochondria_config.dir, train_config.dir];

fprintf('Loading features .. ');
features_mitochondria_trn = [];
features_not_mitochondria_trn = [];
for case_id = train_config.case_ids
  data = load2(sprintf([image_dir, train_config.image_prefix, '.mitochondria_', ...
    feature_config.type, feature_config.version, '.mat'], case_id));

  s1 = size(data.features_mitochondria, 1);
  s2 = size(data.features_not_mitochondria, 1);
  
  index1 = randperm(s1);
  index2 = randperm(s2);
  
  s12 = min(s1,s2);
  
  features_mitochondria_trn = [features_mitochondria_trn; data.features_mitochondria(index1(1:s12),:)];
  features_not_mitochondria_trn = [features_not_mitochondria_trn; data.features_not_mitochondria(index2(1:s12),:)];
end;

if(size(features_mitochondria_trn,1) > MAX_N_VECT)
  skip = max(1, round(features_mitochondria_trn/MAX_N_VECT));
  features_mitochondria_trn = features_mitochondria_trn(1:skip:end, :);
end;
if(size(features_not_mitochondria_trn,1) > MAX_N_VECT)
  skip = max(1, round(features_not_mitochondria_trn/MAX_N_VECT));
  features_not_mitochondria_trn = features_not_mitochondria_trn(1:skip:end, :);
end;
features = [features_mitochondria_trn; features_not_mitochondria_trn];
labels = [ones(size(features_mitochondria_trn,1),1); -ones(size(features_not_mitochondria_trn,1),1)];
fprintf('done\n');


fprintf('Training boosted classifier .. ');
weak_learner = tree_node_w(model_config.tree_depth);
[GLearners, GWeights] = GentleAdaBoost(weak_learner, features', labels', model_config.n_iteration);
window_sizes = data.window_sizes;
intensity_bins = data.intensity_bins;
log = 'histogram equalized; intensity histogram';
MaxIter = model_config.n_iteration;
save2([image_dir, train_config.model_prefix, '.', model_config.type, model_config.version, '.', ...
  feature_config.type, feature_config.version, '.mat'], 'GLearners', 'GWeights', 'MaxIter', ...
  'window_sizes', 'intensity_bins', 'log');
fprintf('done\n');

end
