function to_merge_sets = get_merge_decision_mito_segment(label_map, image, ...
  boundary_sets, merge_criterion_param, image_prefix, config) %#ok<INUSL>
% to_merge_sets = get_merge_decision_region_max(label_map, boundary_map, ...
%   boundary_sets, merge_criterion_param)
%
% Merge criterion for agglomerative prune-classify segmentation. Merges if
% a patch along the boundary between two segments has no value greater than
% a threshold
%
% Shaul Druckmann
% Janelia Farm Research Campus, HHMI.
%

mito = load_mitochondria_detection(...
  [get_reconstruction_dir(config), config.mitochondria.dir, ...
  image_prefix, config.mitochondria.apply.save_suffix, '.mat']);

mito_pixel_threshold = 2;
mito_mask = mito.mitochondria_detection_confidence > mito_pixel_threshold;

mito_perc_threshold = 0.95;

segment_mito_props = regionprops(label_map, mito_mask, 'MeanIntensity');

is_segment_mito = [segment_mito_props(:).MeanIntensity]>mito_perc_threshold;

to_merge_sets = [];
for i = 1:length(boundary_sets)
  if(is_segment_mito(boundary_sets(i).segment_label_0) && ...
      is_segment_mito(boundary_sets(i).segment_label_1))
    to_merge_sets(end+1,:) = ...
      [boundary_sets(i).segment_label_0, boundary_sets(i).segment_label_1]; %#ok<AGROW>
  end
end

return
end
