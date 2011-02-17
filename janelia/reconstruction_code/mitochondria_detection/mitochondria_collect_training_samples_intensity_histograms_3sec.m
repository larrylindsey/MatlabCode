function mitochondria_collect_training_samples_intensity_histograms_3sec(config)
% mitochondria_collect_training_samples_intensity_histograms_3sec(config)
% collect training samples given an image and mitochondria mask
% Features: intensity histograms from k, k-1 and k+1 sections to provide 3D
% context.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  ~Feb.2008  init code
% v1  04112008  modified for reconstruction pipeline
% v2  06062008  adding features from 3 sections - might work better in
%               KHLateral
%

fprintf('Constructing feature vectors for training ..\n');

is_debugging = config.DEBUG;

mitochondria_config = config.mitochondria;
train_config = mitochondria_config.train;
feature_config = mitochondria_config.feature;

version_names{1} = 'heq_intensity_hist_3sec';
version_names{2} = 'intensity_hist_3sec';
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
  img_0 = im2double( imread( sprintf([image_dir, train_config.image_prefix, ...
    train_config.image_suffix], case_id) ) );
  img_p = im2double( imread( sprintf([image_dir, train_config.image_prefix, ...
    train_config.image_suffix], case_id-1) ) );
  img_n = im2double( imread( sprintf([image_dir, train_config.image_prefix, ...
    train_config.image_suffix], case_id+1) ) );
  if(version_id==1)
    img_0 = histeq(img_0);
    img_p = histeq(img_p);
    img_n = histeq(img_n);
  end;
  figure(1); imshow(img_0); hold on;
  title('Training Image');
  
  load2( sprintf([image_dir, train_config.image_prefix, train_config.annot_suffix, '.mat'], ...
    case_id));
  zone_mask = zone_mask>0;
  figure(2); imshow(zone_mask);
  title('Mask of Annotated Mitochondria');
  
  mitochondria_mask = imerode(zone_mask, strel('disk', 1));
  figure(3); imshow(mitochondria_mask);
  title('Mask of Positive Samples');
  
  [features_mitochondria_0, junk] = compute_neighbor_intensity_histogram(img_0, double(mitochondria_mask), ...
    window_sizes, intensity_bins);
  [features_mitochondria_p, junk] = compute_neighbor_intensity_histogram(img_p, double(mitochondria_mask), ...
    window_sizes, intensity_bins);
  [features_mitochondria_n, junk] = compute_neighbor_intensity_histogram(img_n, double(mitochondria_mask), ...
    window_sizes, intensity_bins);
  features_mitochondria = [features_mitochondria_0, features_mitochondria_p, ...
    features_mitochondria_n];
  
  not_mitochondria_mask = imerode(1-zone_mask, strel('disk', 1));
  not_mitochondria_mask_e = imerode(not_mitochondria_mask, strel('disk', 15));
  not_mitochondria_band = not_mitochondria_mask.*(1-not_mitochondria_mask_e);
  
  rand_sample_mask = rand(size(img_0))>0.85;
  
  not_mitochondria_mask = max(not_mitochondria_band, not_mitochondria_mask_e .* rand_sample_mask);
  figure(4); imshow(not_mitochondria_mask);
  title('Mask of Negative Samples');
  
  [features_not_mitochondria_0, junk] = compute_neighbor_intensity_histogram(img_0, double(not_mitochondria_mask), ...
    window_sizes, intensity_bins);
  [features_not_mitochondria_p, junk] = compute_neighbor_intensity_histogram(img_p, double(not_mitochondria_mask), ...
    window_sizes, intensity_bins);
  [features_not_mitochondria_n, junk] = compute_neighbor_intensity_histogram(img_n, double(not_mitochondria_mask), ...
    window_sizes, intensity_bins);
  features_not_mitochondria = [features_not_mitochondria_0, features_not_mitochondria_p, ...
    features_not_mitochondria_n];
  
  save2(sprintf([image_dir, train_config.image_prefix, '.mitochondria_', feature_config.type, ...
    feature_config.version, '.mat'], case_id), 'features_mitochondria', 'features_not_mitochondria', ...
    'window_sizes', 'intensity_bins');
  fprintf('done\n');
end;



