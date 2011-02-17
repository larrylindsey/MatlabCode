function linkage_train_overlap_hist_boost_multi_tile(config)
% linkage_train_overlap_hist_boost_multi_tile(config)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%

if(~isfield(config.linkage, 'is_verbose'))
  config.linkage.is_verbose = true;
end
if(~isfield(config.linkage, 'is_verbose_figures'))
  config.linkage.is_verbose_figures = false;
end
if(config.linkage.is_verbose)
  fprintf('START: linkage_train_overlap_hist_boost_multi_tile\n');
end

stack_config = config.stack;
linkage_config = config.linkage;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
train_config = linkage_config.train;

if(strcmp(model_config.type, 'boost')==0)
  error('Classifier type does not match with called function. Exiting');
end;

if(~isfield(feature_config, 'suffix'))
  feature_config.suffix = '';
end

linkage_dir = get_linkage_dir(config);

if(config.linkage.is_verbose)
  fprintf('Putting together all training feature vectors of all sections\n');
end
case_id = stack_config.case_ids(1);
fprintf('case_id: %d\n', case_id);
image_prefixes_1 = get_image_prefixes_subdirs(config, case_id);

features = [];
link_labels = [];
for i = 2:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  fprintf('case_id: %d\n', case_id);

  image_prefixes_2 = get_image_prefixes_subdirs(config, case_id);
  
  % compute a linkage graph for each pair of overlapping tiles
  for tile_1 = 1:length(image_prefixes_1)
    for tile_2 = 1:length(image_prefixes_2)
      fprintf('--- Tile pair ---\n%d,%d\n', tile_1, tile_2);
      fprintf('%s\n%s\n', image_prefixes_1{tile_1}, image_prefixes_2{tile_2});
      
      file_name_prefix = get_file_name_from_tuple(...
        [get_reconstruction_dir(config), linkage_config.dir], ...
        image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'flp.');
      file_name_suffix = [feature_config.type, feature_config.suffix];
      file_name_features = [file_name_prefix, file_name_suffix, ...
        config.segmentation_choose.choice.seg_suffix, '.mat'];
      try
        if(linkage_config.is_verbose)
          fprintf('Attempting to load ground-truth:\n%s\n', file_name_features);
        end
        linkage_features = load2(file_name_features, 'features', ...
          'link_labels');
      catch %#ok<CTCH>
        if(linkage_config.is_verbose)
          fprintf('failed, skipping.\n');
        end
        continue;
      end
      
      features = [features; linkage_features.features]; %#ok<AGROW>
      link_labels = [link_labels; linkage_features.link_labels]; %#ok<AGROW>
    end
  end

  image_prefixes_1 = image_prefixes_2;
end;

if(isempty(features))
  error('There were no training samples. Weird!\n');
end

fprintf('Collected %d feature vectors of dimensionality %d\n', ...
  size(features,1), size(features,2));

save2([linkage_dir, 'link_train_feat', '.', ...
  feature_config.type, feature_config.suffix, train_config.save_suffix, '.mat'], ...
  'features', 'link_labels', 'feature_config');

%
% train a adaboost based classifier
%
relative_error_weights = ones(size(link_labels));
if(isfield(model_config, 'merge_error_relative_weight'))
  relative_error_weights(link_labels<0) = ...
    model_config.merge_error_relative_weight;
end

fprintf('Training classifier for linkage across sections ..\n');
MaxIter = model_config.n_iteration;
tree_depth = model_config.tree_depth; % 25,2
weak_learner = tree_node_w(tree_depth);
[classifier.GLearners, classifier.GWeights] = ...
  GentleAdaBoost(weak_learner, features', link_labels', MaxIter, ...
  relative_error_weights');
fprintf('done.\n');

n_dimension = size(features,2); %#ok<NASGU>
save2([linkage_dir, 'link_model', '.', ...
  model_config.type, model_config.suffix, '.', ...
  feature_config.type, feature_config.suffix, ...
  train_config.save_suffix, '.mat'], ...
  'classifier', 'n_dimension', 'model_config', 'feature_config');
fprintf('done\n');


if(linkage_config.is_verbose)
  fprintf('STOP: linkage_train_overlap_hist_boost_multi_tile\n');
end
return
end
