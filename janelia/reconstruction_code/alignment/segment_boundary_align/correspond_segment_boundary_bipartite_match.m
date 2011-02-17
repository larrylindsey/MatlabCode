function [segment_correspondence, votes_border] = ...
  correspond_segment_boundary_bipartite_match...
  (segment_0, border_mask_0, stitch_segments_0, segment_1, border_mask_1, config)
% segment_correspondence = correspond_segment_boundary_bipartite_match(segment_0, ...
%    stitch_segments_0, segment_1, config)
% Correspond segments from two overlapping segment maps. The map in
% segment_0 is assumed to lay "above" the one from segment_1. The stitching
% is performed for the set of segments in segment_0 specified in
% stitch_segments_0.
% Algorithm: match the boundary elements of segments using bipartite
% matching (CSA++ Berkeley segmentation evaluation). Each match is
% considered a vote of correspondence between two segments
%
% Input:
%   segment_0, segment_1    Segment maps with integer labels for segments,
%                           0's for boundaries and -1's for padded pixels.
%   border_mask_0,1         Masks of borders to be considered for
%                             correspondence
%   stitch_segments_0       Set of segments in segment_0 that should be
%                           stitched with segments in segment_1.
%   config                  Reconstruction config struct
% Output:
%   segment_correspondence  pairing of correspondences (2xN matrix)
%

align_config = config.align_segmentation;
if(~isfield(align_config, 'is_verbose_figures'))
  align_config.is_verbose_figures = false;
end
combined_mask = border_mask_0 .* border_mask_1;
boundary_map_0 = double(combined_mask .* (segment_0==0) .* ...
  (imdilate(segment_0, strel('disk',2))>0));
boundary_map_1 = double(combined_mask .* (segment_1==0) .* ...
  (imdilate(segment_1, strel('disk',2))>0));

if(max(boundary_map_0(:))<=0 || max(boundary_map_1(:))<=0)
  segment_correspondence = [];
  return;
end

delta_margin = align_config.delta_margin/sqrt(sum(size(boundary_map_0).^2));

[match_0, match_1] = correspondPixels(boundary_map_0, ...
  boundary_map_1, delta_margin); %#ok<NASGU>
if(align_config.is_verbose_figures)
  figure;
  di = repmat(boundary_map_0, [1 1 3]);
  di(:,:,3) = boundary_map_1;
  imshow(di);
  [match_y_0, match_x_0] = find(match_0);
  [match_y_1, match_x_1] = ind2sub(size(boundary_map_0), match_0(match_0>0));
  hold on;
  plot([match_x_0'; match_x_1'], [match_y_0'; match_y_1'], 'r');
  hold off;
end

b_to_b_map = [find(match_0>0)'; match_0(match_0>0)'];
votes = correspond_segments_from_boundary_correspondence(b_to_b_map, ...
  segment_0, segment_1);
% votes_border = votes(:, ismember(votes(1,:), stitch_segments_0));
votes_border = votes(:, votes(1,:)>0 & votes(2,:)>0);

% Apply minimum threshold
votes_significant = votes_border(:, votes_border(3,:)>=align_config.min_n_votes);

% Of the significant correspondences, each segment of each segment-map is
% made to have atleast mapping with highest number of votes. This will
% enable mapping.
[sorted_n_votes, sort_id] = sort(votes_significant(3,:)); %#ok<ASGLU>
votes_significant_sorted = votes_significant(:,sort_id);

[label_0, back_ref_0] = unique(votes_significant_sorted(1,:), 'last'); %#ok<ASGLU>
segment_correspondence_0 = votes_significant_sorted(:, back_ref_0);

[label_1, back_ref_1] = unique(votes_significant_sorted(2,:), 'last'); %#ok<ASGLU>
segment_correspondence_1 = votes_significant_sorted(:, back_ref_1);

segment_correspondence = ...
  unique([segment_correspondence_0'; segment_correspondence_1'], 'rows')';

return
end
