function mitochondria_apply_boosted_detector_intensity_norm_hist_3s2(config)
% mitochondria_apply_boosted_detector_intensity_norm_hist_3s2(config)
% apply the trained mitochondria detector to confidence maps on image plane
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  ~Feb.2008  init code
% v1  04112008  modified for reconstruction pipeline
%

stack_config = config.stack;
mitochondria_config = config.mitochondria;
train_config = mitochondria_config.train;
feature_config = mitochondria_config.feature;
model_config = mitochondria_config.model;
apply_config = mitochondria_config.apply;

if(strcmp(model_config.type, 'boost')==0)
  error('Classifier type does not match with called function. Exiting');
end;
if(~isfield(model_config, 'version'))
  model_config.version = '';
end

version_names{1} = 'heq_intensity_norm_hist_3s2';
version_names{2} = 'intensity_norm_hist_3s2';
version_id = find(strcmp(feature_config.type, version_names));
if(isempty(version_id))
  error('Feature type does not match with called function. Exiting');
end;
if(~isfield(feature_config, 'version'))
  feature_config.version = '';
end

load2([get_reconstruction_dir(config), mitochondria_config.dir, train_config.dir, ...
  train_config.model_prefix, '.', model_config.type, model_config.version, '.', ...
  feature_config.type, feature_config.version, '.mat'], 'GLearners', 'GWeights', 'MaxIter', ...
  'window_sizes', 'intensity_bins');

for case_id = stack_config.case_ids
  
  fprintf('case %d: ', case_id);
  fprintf('collecting samples..');
  [img_0, ip, isb, is_to_be_processed] = ...
    get_image_from_stack(config, case_id);
  if(is_to_be_processed(1)==0)
    fprintf('is_to_be_processed=false, skipping\n');
    continue;
  end
    
  img_0 = img_0{1};
  img_p = get_image_from_stack(config, case_id-1);
  img_p = img_p{1};
  img_n = get_image_from_stack(config, case_id+1);
  img_n = img_n{1};
  
  if(~isempty(stack_config.roi))
    img_0 = img_0(stack_config.roi.ymin:stack_config.roi.ymax, ...
      stack_config.roi.xmin:stack_config.roi.xmax);
    img_p = img_p(stack_config.roi.ymin:stack_config.roi.ymax, ...
      stack_config.roi.xmin:stack_config.roi.xmax);
    img_n = img_n(stack_config.roi.ymin:stack_config.roi.ymax, ...
      stack_config.roi.xmin:stack_config.roi.xmax);
  end;
  
  if(version_id==1)
    img_0 = histeq(img_0);
    img_p = histeq(img_p);
    img_n = histeq(img_n);
  end;
  
  temp = [];
  temp(:,:,1) = img_0;
  temp(:,:,2) = img_p;
  temp(:,:,3) = img_n;
  
  img_1 = min(temp, [], 3);
  img_2 = mean(temp, 3);
  img_3 = max(temp, [], 3);
  
  sample_mask = ones(size(img_0));
  [test_features_1, sample_loc] = compute_neighbor_intensity_histogram_norm(img_1, double(sample_mask), ...
    window_sizes, intensity_bins);
  [test_features_2, sample_loc] = compute_neighbor_intensity_histogram_norm(img_2, double(sample_mask), ...
    window_sizes, intensity_bins);
  [test_features_3, sample_loc] = compute_neighbor_intensity_histogram_norm(img_3, double(sample_mask), ...
    window_sizes, intensity_bins);
  test_features = [test_features_1, test_features_2, test_features_3];
  fprintf('done. ');

  fprintf('classifying samples..');
  weak_learner = tree_node_w(model_config.tree_depth);
  d = Classify(GLearners, GWeights, test_features');
  fprintf('done. ');

  mitochondria_detection_confidence = zeros(size(img_0))-inf;
  sample_loc_id = (sample_loc(:,1)-1)*size(img_0,1) + sample_loc(:,2);
  mitochondria_detection_confidence(sample_loc_id) = d;
  
  fprintf('saving ..');
  save2(sprintf([get_reconstruction_dir(config), mitochondria_config.dir, stack_config.image_prefix, ...
    apply_config.save_suffix, '.mat'], case_id), 'mitochondria_detection_confidence');
  fprintf('done. ');
  
  fprintf('\n');
end;
