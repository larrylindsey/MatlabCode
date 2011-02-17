function body_to_merge = get_linkage_mergers(z_i, z_j, body_id_i, body_ids_j)
% body_to_merge = get_linkage_mergers(z1, z2, body_id_i, body_ids_j)
% Compute the bodies among body_ids_j that should be merged with body_id_i
% through 3D linkage confidence.
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

global reconstruction_config

switch(reconstruction_config.linkage.feature.type)
  case {'intensity_pair_hist'}
    switch(reconstruction_config.linkage.model.type)
      case 'boost'
        body_to_merge = get_linkage_mergers_intensity_pair_boost(z_i, z_j, body_id_i, body_ids_j);
      otherwise
        error('Linkage model type not recognized');
    end
  case {'intensity_pair_hist_v2', 'intensity_pair_hist_v2b', ...
      'intensity_pair_hist_v2c', 'intensity_pair_hist_v2d'}
    switch(reconstruction_config.linkage.model.type)
      case 'boost'
        body_to_merge = get_linkage_mergers_intensity_pair_v2_boost(z_i, z_j, body_id_i, body_ids_j);
      otherwise
        error('Linkage model type not recognized');
    end
  otherwise
    error('Linkage feature type not recognized');
end

return;
end
