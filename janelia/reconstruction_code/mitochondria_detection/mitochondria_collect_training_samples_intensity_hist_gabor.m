function mitochondria_collect_training_samples_intensity_hist_gabor(config)
% mitochondria_collect_training_samples_intensity_hist_gabor(config)
% collect training samples given an image and mitochondria mask
%  - intensity and Gabor-filter histograms within concentric patches
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  ~Feb.2008  init code
% v1  04112008  modified for reconstruction pipeline
% v2  05012008  borrowed from
% mitochondria_collect_training_samples_intensity_histograms.m added Gabor filter 
%

fprintf('Constructing feature vectors for training ..\n');

mitochondria_config = config.mitochondria;
train_config = mitochondria_config.train;
feature_config = mitochondria_config.feature;

if(isfield(mitochondria_config, 'verbose'))
  is_verbose = mitochondria_config.verbose;
else
  is_verbose = true;
end;

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

image_dir = [get_reconstruction_dir(config), mitochondria_config.dir, train_config.dir];

window_sizes = feature_config.window_sizes';
intensity_bins = feature_config.intensity_bins';

for case_id = train_config.case_ids
  fprintf('Collecting samples for image %d ..', case_id);
  img = im2double( imread( sprintf([image_dir, train_config.image_prefix, ...
    train_config.image_suffix], case_id) ) );
  img = img/max(img(:));
  img_gabor = img;
  if(version_id==1)
    img = histeq(img);
  end;
  if(version_id==3)
    img = histeq(img);
    img_gabor = img;
  end;
  if(is_verbose)
    figure(1); clf; imshow(img);
    title('Training Image');
  end;
  
  load2( sprintf([image_dir, train_config.image_prefix, train_config.annot_suffix, '.mat'], ...
    case_id));
  zone_mask = zone_mask>0;
  if(is_verbose)
    figure(2); imshow(zone_mask);
    title('Mask of Annotated Mitochondria');
  end
  
  mitochondria_mask = imerode(zone_mask, strel('disk', 1));
  if(is_verbose)
    figure(3); imshow(mitochondria_mask);
    title('Mask of Positive Samples');
  end
  
  not_mitochondria_mask = imerode(1-zone_mask, strel('disk', 1));
  not_mitochondria_mask_e = imerode(not_mitochondria_mask, strel('disk', 15));
  not_mitochondria_band = not_mitochondria_mask.*(1-not_mitochondria_mask_e);
  
  rand_sample_mask = rand(size(img))>0.85;
  
  not_mitochondria_mask = max(not_mitochondria_band, not_mitochondria_mask_e .* rand_sample_mask);
  if(is_verbose)
    figure(4); imshow(not_mitochondria_mask);
    title('Mask of Negative Samples');
  end

  % Intensity hist
  [features_mitochondria, junk] = compute_neighbor_intensity_histogram(img, double(mitochondria_mask), ...
    window_sizes, intensity_bins);
  
  [features_not_mitochondria, junk] = compute_neighbor_intensity_histogram(img, double(not_mitochondria_mask), ...
    window_sizes, intensity_bins);
  
  % Gabor filter
  for f = feature_config.gabor.freqs
    map = 1-img_gabor;
    
    [g1,g2] = gaborfilter_maxTheta(map, feature_config.gabor.size, feature_config.gabor.size, ...
      f, feature_config.gabor.n_theta);
    
    max_g1 = max(g1(:)); min_g1 = min(g1(:)); g1 = (g1-min_g1)/(max_g1-min_g1+eps);
    [feat, junk] = compute_neighbor_intensity_histogram(g1, double(mitochondria_mask), ...
      window_sizes, intensity_bins);
    features_mitochondria = [features_mitochondria, feat];
    [feat, junk] = compute_neighbor_intensity_histogram(g1, double(not_mitochondria_mask), ...
      window_sizes, intensity_bins);
    features_not_mitochondria = [features_not_mitochondria, feat];
    
    max_g2 = max(g2(:)); min_g2 = min(g2(:)); g2 = (g2-min_g2)/(max_g2-min_g2+eps);
    [feat, junk] = compute_neighbor_intensity_histogram(g2, double(mitochondria_mask), ...
      window_sizes, intensity_bins);
    features_mitochondria = [features_mitochondria, feat];
    [feat, junk] = compute_neighbor_intensity_histogram(g2, double(not_mitochondria_mask), ...
      window_sizes, intensity_bins);
    features_not_mitochondria = [features_not_mitochondria, feat];
  end
  
  save2(sprintf([image_dir, train_config.image_prefix, '.mitochondria_', feature_config.type, ...
    feature_config.version, '.mat'], case_id), 'features_mitochondria', 'features_not_mitochondria', ...
    'window_sizes', 'intensity_bins');
  fprintf('done\n');
end;
