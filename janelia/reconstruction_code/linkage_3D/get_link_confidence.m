function link_confidence = get_link_confidence(z1, z2, body_id1, body_id2)
% link_confidence = get_link_confidence(z1, z2, body_id1, body_id2)
% Compute the 3D linkage confidence between two bodies, namely body_id_1 and
% body_id_2.
%
% This is to be used for updating the 3D links between bodies during
% proofreading. The function is implemented to be called from gui.m. Makes
% assumptions about existence of al, cat, superpixel_2_seg_map, proof
% variables. 
% The bodies are assumed to be cat->superpixel_2_seg_map->proof.pmap. (See
% gui.m for details).
% Parameters for the linkage are obtained from config.linkage. See
% reconstruction_ex_Alex_volume.m for illustration.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  05282008  init code
%

global al cat superpixel_2_seg_map proof reconstruction_config linkage_model

if(exist('reconstruction_config', 'var')==0 || ~isfield(reconstruction_config, 'linkage'))
  link_confidence = 0
  return;
end

stack_config = reconstruction_config.stack;
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

% if(exist('linkage_model', 'var')==0)
  linkage_model = load([get_reconstruction_dir(reconstruction_config), linkage_config.dir, 'link_model', '.', ...
    model_config.type, '.', feature_config.type, train_config.save_suffix, '.mat'], ...
    'GLearners', 'GWeights', 'MaxIter', 'tree_depth', 'annotation_file', 'intensity_bins');
  weak_learner = tree_node_w(linkage_model.tree_depth);
% end

% get the overlap area
map1 = proof.pmap(1+superpixel_2_seg_map{z1}(1+cat{z1})) == body_id1;
map2 = proof.pmap(1+superpixel_2_seg_map{z2}(1+cat{z2})) == body_id2;
overlap_map = map1 & map2;

% compute the features
image_1 = histeq(im2double(al{z1}));
image_2 = histeq(im2double(al{z2}));
intensity_pair_set = double([image_1(overlap_map(:)), image_2(overlap_map(:))]);

h = hist2(intensity_pair_set, linkage_model.intensity_bins, linkage_model.intensity_bins);

% if the model has been trained of images of resolution different from
% the current one, then the area and length should be adjusted. This
% would be especially useful during bootstrapping.
% h = h * stack_config.area_factor;

intensity_pair_hist_feature = reshape(h, [1, numel(h)]);

features = intensity_pair_hist_feature;

% apply classifier and return
link_confidence = Classify(linkage_model.GLearners, linkage_model.GWeights, features')

return;

end