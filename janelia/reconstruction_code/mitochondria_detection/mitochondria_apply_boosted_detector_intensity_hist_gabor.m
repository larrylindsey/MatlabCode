function mitochondria_apply_boosted_detector_intensity_hist_gabor(config)
% mitochondria_apply_boosted_detector_intensity_hist_gabor(config)
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

weak_learner = tree_node_w(model_config.tree_depth);

if(strcmp(model_config.type, 'boost')==0)
  error('Classifier type does not match with called function. Exiting');
end;
if(~isfield(model_config, 'version'))
  model_config.version = '';
end

version_names{1} = 'heq_intensity_hist_gabor';
version_names{2} = 'intensity_hist_gabor';
version_names{3} = 'heq_intensity_hist_heq_gabor';
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
  [images, image_prefixes, image_sub_dirs, is_to_be_processed] = ...
    get_image_from_stack(config, case_id);
  
  for tile_id = 1:length(images)
    fprintf('tile %d ', tile_id);
    img = images{tile_id};
    image_prefix = image_prefixes{tile_id};
    image_sub_dir = image_sub_dirs{tile_id};
  
    if(is_to_be_processed(tile_id)==0)
      fprintf('is_to_be_processed=false, skipping\n');
      continue;
    end
    
    if(~isempty(stack_config.roi))
      img = img(stack_config.roi.ymin:stack_config.roi.ymax, ...
        stack_config.roi.xmin:stack_config.roi.xmax);
    end;

    img_gabor = img;
    if(version_id==1)
      img = histeq(img);
    end;
    if(version_id==3)
      img = histeq(img);
      img_gabor = img;
    end;
    mitochondria_detection_confidence = zeros(size(img))-inf;

    fprintf('Computing Gabor responses ..');
    g1 = [];
    g2 = [];
    for f = 1:length(feature_config.gabor.freqs)
      fprintf('%d ', f);
      [g1(:,:,f),g2(:,:,f)] = gaborfilter_maxTheta(1-img_gabor, feature_config.gabor.size, feature_config.gabor.size, ...
        feature_config.gabor.freqs(f), feature_config.gabor.n_theta);

      max_g1 = max(max(squeeze(g1(:,:,f)))); min_g1 = min(min(squeeze(g1(:,:,f))));
      g1(:,:,f) = (g1(:,:,f)-min_g1)/(max_g1-min_g1+eps);

      max_g2 = max(max(squeeze(g2(:,:,f)))); min_g2 = min(min(squeeze(g2(:,:,f))));
      g2(:,:,f) = (g2(:,:,f)-min_g2)/(max_g2-min_g2+eps);
    end;
    fprintf('done. ');

    fprintf('Computing confidence values ..');
    band_t = 1:200:size(img,1);
    band_b = min(band_t+199, size(img,1));
    for band = 1:length(band_t);
      sample_mask = zeros(size(img));
      sample_mask(band_t(band):band_b(band), :) = 1;

      % Intensity hist
      [test_features, sample_loc] = compute_neighbor_intensity_histogram(img, double(sample_mask), window_sizes, intensity_bins);

      if(isempty(sample_loc))
        continue;
      end;

      % Gabor filter
      for f = 1:length(feature_config.gabor.freqs)
        [feat, junk] = compute_neighbor_intensity_histogram(g1(:,:,f), double(sample_mask), ...
          window_sizes, intensity_bins);
        test_features = [test_features, feat];

        [feat, junk] = compute_neighbor_intensity_histogram(g2(:,:,f), double(sample_mask), ...
          window_sizes, intensity_bins);
        test_features = [test_features, feat];
      end

      d = Classify(GLearners, GWeights, test_features');

      sample_loc_id = (sample_loc(:,1)-1)*size(img,1) + sample_loc(:,2);
      mitochondria_detection_confidence(sample_loc_id) = d;
    end
    fprintf('done. ');

    fprintf('saving .. ');
    check_for_dir([get_reconstruction_dir(config), mitochondria_config.dir, ...
      image_sub_dir]);
    save_mitochondria_detection(...
      [get_reconstruction_dir(config), mitochondria_config.dir, image_prefix, ...
      apply_config.save_suffix, '.mat'], mitochondria_detection_confidence);
    fprintf('done. ');

    fprintf('\n');
  end
end;

return
end