function linkage_preprocess_image()
% linkage_preprocess_image()
% Preprocess images for computing linkage features during proofreading.
% To be used with get_linkage_mergers().
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  05292008  init code
%

fprintf('Preprocessing images for linkage ... ');
global al reconstruction_config preprocessed_images_linkage

if(exist('reconstruction_config', 'var')==0 || ~isfield(reconstruction_config, 'linkage'))
  preprocessed_images_linkage = [];
  return;
end

linkage_config = reconstruction_config.linkage;
feature_config = linkage_config.feature;

version_names{1} = 'intensity_pair_hist';
version_names{2} = 'intensity_pair_hist_v2';
version_names{3} = 'intensity_pair_hist_v2b';
version_names{4} = 'intensity_pair_hist_v2c';
version_names{5} = 'intensity_pair_hist_v2d';
version_id = find(strcmp(feature_config.type, version_names));
if(isempty(version_id))
  error('Feature type does not match with called function. Exiting');
end;

preprocessed_images_linkage = [];
for i = 1:length(al)
  preprocessed_images_linkage{i} = histeq(im2double(al{i}));
end
fprintf('done\n');

return
end
