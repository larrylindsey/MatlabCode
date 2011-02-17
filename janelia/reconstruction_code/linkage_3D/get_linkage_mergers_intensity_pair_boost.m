function body_to_merge = get_linkage_mergers_intensity_pair_boost(z_i, z_j, body_id_i, body_ids_j)
% body_to_merge = get_linkage_mergers_intensity_pair_boost(z1, z2, body_id_i, body_ids_j)
% Compute the bodies among body_ids_j that should be merged with body_id_i
% through 3D linkage confidence.
% Features: pairs of pixel intensities from overlapping segment regions.
%
% This is to be used for updating the 3D links between bodies during
% proofreading. The function is implemented to be called from gui.m. Makes
% assumptions about existence of preprocessed_images_linkage, cat,
% superpixel_2_seg_map, proof variables.
%
% preprocessed_images_linkage holds preprocessed images for faster feature
% computation. See linkage_preprocess_image.m for details.
%
% The bodies are assumed to be cat->superpixel_2_seg_map->proof.pmap. (See
% gui.m for details).
%
% Parameters for the linkage are obtained from config.linkage. See
% reconstruction_ex_Alex_volume.m for illustration.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  05282008  init code
%

global cat superpixel_2_seg_map proof
global reconstruction_config preprocessed_images_linkage link_threshold

if(exist('reconstruction_config', 'var')==0 || ~isfield(reconstruction_config, 'linkage'))
  link_confidence = 0;
  return;
end

linkage_config = reconstruction_config.linkage;
train_config = linkage_config.train;
feature_config = linkage_config.feature;
model_config = linkage_config.model;

if(strcmp(feature_config.type, 'intensity_pair_hist')==0)
  error('Feature type does not match with called function. Exiting');
end;
if(strcmp(model_config.type, 'boost')==0)
  error('Classifier type does not match with called function. Exiting');
end;

if(exist('linkage_model', 'var')==0)
  global linkage_model
  if(isfield(reconstruction_config, 'link_model_file'))
    linkage_model = load(reconstruction_config.link_model_file, ...
    'GLearners', 'GWeights', 'MaxIter', 'tree_depth', 'annotation_file', 'intensity_bins');
  else
    linkage_model = load([get_reconstruction_dir(reconstruction_config), linkage_config.dir, 'link_model', '.', ...
      model_config.type, '.', feature_config.type, train_config.save_suffix, '.mat'], ...
      'GLearners', 'GWeights', 'MaxIter', 'tree_depth', 'annotation_file', 'intensity_bins');
  end
  weak_learner = tree_node_w(linkage_model.tree_depth);
end

% compute the features from the overlap area
map_i = proof.pmap(1+superpixel_2_seg_map{z_i}(1+cat{z_i})) == body_id_i;
body_map = proof.pmap(1+superpixel_2_seg_map{z_j}(1+cat{z_j}));
body_to_merge = [];
for b_j = body_ids_j
  map_j = body_map == b_j;
  overlap_map = map_i & map_j;

  intensity_pair_set = double([preprocessed_images_linkage{z_i}(overlap_map(:)), ...
    preprocessed_images_linkage{z_j}(overlap_map(:))]);

  h = hist2(intensity_pair_set, linkage_model.intensity_bins, linkage_model.intensity_bins);

  % if the model has been trained of images of resolution different from
  % the current one, then the area and length should be adjusted. This
  % would be especially useful during bootstrapping.
  % h = h * stack_config.area_factor;

  intensity_pair_hist_feature = reshape(h, [1, numel(h)]);

  features = intensity_pair_hist_feature;

  % apply classifier and return
  link_confidence = Classify(linkage_model.GLearners, linkage_model.GWeights, features');
  
  if(link_confidence>link_threshold)
    body_to_merge(end+1) = b_j;
  end;
end

return;
end
