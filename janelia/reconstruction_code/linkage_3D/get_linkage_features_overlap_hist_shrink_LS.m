function [features, label_pairs, label_pairs_original] = ...
  get_linkage_features_overlap_hist_shrink_LS(...
  label_map_1, image_1, label_map_2, image_2, transforms_tp, transforms_tp_rev, ...
  feature_config, label_map_original_1, label_map_original_2)
% get_linkage_features_overlap_hist_shrink_LS(config)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

if(nargin<8)
  label_map_original_1 = [];
  label_map_original_2 = [];
end

region_overlap_features = [];
region_overlap_features_rev = [];
region_overlap_features_nominal = [];
region_overlap_features_rev_nominal = [];
if(~( isempty(label_map_1) || ...
    isempty(label_map_2) || ...
    isempty(transforms_tp) || ...
    isempty(transforms_tp.transforms) || ...
    isempty(transforms_tp.map_mask) || ...
    size(transforms_tp.transforms, 2)<=4 || ...
    size(transforms_tp_rev.transforms, 2)<=4))
  
  fprintf('Computing linkage features between segment maps ...\n');
  region_overlap_features = ...
    collect_seg_overlap_stats_interior_hist_paffine(...
    uint32(label_map_1), ...
    uint32(label_map_2), ...
    uint8(255*image_1), uint8(255*image_2), ...
    uint8(feature_config.bin_size), ...
    transforms_tp.map_mask, transforms_tp.transforms);
  if(~isempty(label_map_original_1))
    region_overlap_features_nominal = ...
      collect_seg_overlap_stats_interior_hist_paffine(...
      uint32(label_map_original_1), ...
      uint32(label_map_original_2), ...
      uint8(255*image_1), uint8(255*image_2), ...
      uint8(feature_config.bin_size), ...
      transforms_tp.map_mask, transforms_tp.transforms);
  else
    region_overlap_features_nominal = [];
  end
end
if(~( isempty(label_map_1) || ...
    isempty(label_map_2) || ...
    isempty(transforms_tp_rev) || ...
    isempty(transforms_tp_rev.transforms) || ...
    isempty(transforms_tp_rev.map_mask) || ...
    size(transforms_tp_rev.transforms, 2)<=4 || ...
    size(transforms_tp.transforms, 2)<=4))
  
  fprintf('Computing linkage features between segment maps in reverse...\n');
  region_overlap_features_rev = ...
    collect_seg_overlap_stats_interior_hist_paffine(...
    uint32(label_map_2), ...
    uint32(label_map_1), ...
    uint8(255*image_2), uint8(255*image_1), ...
    uint8(feature_config.bin_size), ...
    transforms_tp_rev.map_mask, transforms_tp_rev.transforms);
  if(~isempty(label_map_original_1))
    region_overlap_features_rev_nominal = ...
      collect_seg_overlap_stats_interior_hist_paffine(...
      uint32(label_map_original_2), ...
      uint32(label_map_original_1), ...
      uint8(255*image_2), uint8(255*image_1), ...
      uint8(feature_config.bin_size), ...
      transforms_tp_rev.map_mask, transforms_tp_rev.transforms);
  else
    region_overlap_features_rev_nominal = [];
  end
end
      
region_features_1 = [];
region_features_2 = [];
if(~isempty(region_overlap_features) || ...
    ~isempty(region_overlap_features_rev))
  interior_hist_t = collect_segment_stats_interior_hist(...
    uint32(label_map_1), ...
    uint8(255*image_1), uint8(feature_config.bin_size));
  interior_hist = [];
  interior_hist(interior_hist_t(:,1)+1, :) = interior_hist_t(:, 2:end);
  region_features_1 = interval_sum(interior_hist);
  
  interior_hist_t = collect_segment_stats_interior_hist(...
    uint32(label_map_2), ...
    uint8(255*image_2), uint8(feature_config.bin_size));
  interior_hist = [];
  interior_hist(interior_hist_t(:,1)+1, :) = interior_hist_t(:, 2:end);
  region_features_2 = interval_sum(interior_hist);
end

if(~isempty(region_overlap_features_rev))
  region_overlap_features = [region_overlap_features; ...
    region_overlap_features_rev(:, [2 1 3:end])];
end

if(isempty(region_overlap_features))
  features = [];
  label_pairs = [];
else
  features = [interval_sum(region_overlap_features(:, 3:end)), ...
    region_features_1(1+region_overlap_features(:,1), :), ...
    region_features_2(1+region_overlap_features(:,2), :)];
  
  label_pairs = region_overlap_features(:,1:2);
end

label_pairs_original = [];
if(~isempty(region_overlap_features_nominal))
  label_pairs_original = region_overlap_features_nominal(:, [1 2]);
end
if(~isempty(region_overlap_features_rev_nominal))
  label_pairs_original = [label_pairs_original; ...
    region_overlap_features_rev_nominal(:, [2 1])];
end

return
end