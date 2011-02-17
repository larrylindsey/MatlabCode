function mitochondria_segment_train_classifier(config)
% mitochondria_segment_boost_train(config)
% Train boosted classifier to give model for computing mitochondria
% confidence values for segments.
%
% Nicholas Sofroniew,
% Univ. of Cambridge, UK.
% Visitor March 2009, JFRC, HHMI.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  03132009  init code
%

stack_config = config.stack;
mitochondria_config = config.mitochondria;
segment_config = mitochondria_config.segment;
classify_config = segment_config.classify;
feature_config = classify_config.feature;
model_config = classify_config.model;

if(~isfield(mitochondria_config, 'is_verbose'))
  mitochondria_config.is_verbose = true;
end

mito_dir = [get_reconstruction_dir(config), mitochondria_config.dir];
if(mitochondria_config.is_verbose)
  fprintf('\nSTART: mitochondria_segment_train_classifier\n');
end

if(mitochondria_config.is_verbose)
  fprintf('Loading segment features and ground truth ... ');
end

features_stack=[];
gt_stack=[];
for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  if(mitochondria_config.is_verbose)
    fprintf('plane:%d\n', case_id);
  end
  image_prefixes = get_image_prefixes_subdirs(config, case_id);
  image_prefix = image_prefixes{1};

  load_file_name = [mito_dir, image_prefix, ...
    '.mito.seg', segment_config.mito_seg_suffix, ...
    '.features', feature_config.features_suffix, ...
    '.mat'];
  segment_features = load2(load_file_name, 'features');
  features_stack=[features_stack; segment_features.features ];

  load_file_name = [mito_dir, image_prefix, ...
    '.mito.seg', segment_config.mito_seg_suffix, ...
    '.gt', segment_config.gt_suffix, ...
    '.mat'];
  ground_truth_labels = load2(load_file_name, 'seg_labels');
  gt_stack=[gt_stack; ground_truth_labels.seg_labels];

end
if(mitochondria_config.is_verbose)
  fprintf('done.\n');
end

%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
%model.tree_depth = 2;
%model.n_iteration = 50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(mitochondria_config.is_verbose)
  fprintf('Training classifier ... ');
end
labels_train=2*gt_stack-1;
samples_train=[features_stack(:,2:end)];
weak_learner=tree_node_w(model_config.tree_depth);

%%%% weight segments according to their area using features_stack(:,1)
[model_config.GLearners, model_config.GWeights]=GentleAdaBoost(weak_learner, samples_train',labels_train',model_config.n_iteration, features_stack(:,1)');
model_config.weak_learner = weak_learner;

if(mitochondria_config.is_verbose)
  fprintf('done.\n');
end

%% Save model...
save_file_name = [mito_dir, ...
  '.mito.seg', segment_config.mito_seg_suffix, ...
  '.features', feature_config.features_suffix, ...
  '.gt', segment_config.gt_suffix, ...
  '.model.', num2str(model_config.tree_depth), ...
  '.',num2str(model_config.n_iteration), ...
  '.mat'];
save2(save_file_name, 'model_config');

if(mitochondria_config.is_verbose)
  fprintf('\nSTOP: mitochondria_segment_train_classifier\n');
end
return;
end
