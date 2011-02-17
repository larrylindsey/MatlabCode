function [n_split, n_merge, split_segment_id, merge_segment_id] = ...
  compute_split_merge_error(seg_2_body_gt, seg_2_body)
% [n_split, n_merge] = compute_split_merge_error(seg_2_body_map_gt, seg_2_body_map)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%

seg_2_body_map_gt = zeros(max([seg_2_body_gt(:,1); seg_2_body(:,1)])+1, 1);
seg_2_body_map_gt(1+seg_2_body_gt(:,1)) = seg_2_body_gt(:,2);

seg_2_body_map = zeros(max([seg_2_body_gt(:,1); seg_2_body(:,1)])+1, 1);
seg_2_body_map(1+seg_2_body(:,1)) = seg_2_body(:,2);

% In the comments and variables, "as" refers to automatic segmentation and
% "gt" refers to ground truth.

% Assign each auto seg. body to one ground-truth body based on maximal
% overlap
[label_pair, count] = count_row_occurence([seg_2_body_map, seg_2_body_map_gt]);
count = count(label_pair(:,1)>0 & label_pair(:,2)>0, :);
label_pair = label_pair(label_pair(:,1)>0 & label_pair(:,2)>0, :);
as_label_overlap = sortrows([count, label_pair], [1 2]);
[junk, I] = unique(as_label_overlap(:, 2), 'rows', 'last');
as_to_gt_map = as_label_overlap(I, [2 3]);
as_to_gt_overlap = as_label_overlap(I,1);

% Number of segments that need to split away from their bodies and assigned
% to some other body is the volume of auto. seg. body - volume of overlap
% of the auto. seg. body with the assigned gt. body.

% Find volumes of auto. seg bodies.
[as_label, count] = count_row_occurence(seg_2_body_map);
count = count(as_label>0);
as_label = as_label(as_label>0);
as_volume = [];
as_volume(as_label+1) = count;

as_gt_overlap = zeros(1, max(as_to_gt_map(:,1))+1);
as_gt_overlap(as_to_gt_map(:,1)+1) = as_to_gt_overlap;

n_split = sum(as_volume - as_gt_overlap);

% List of segments that are not in correct bodies and must be split away
split_segment_id = ...
  find(seg_2_body_map_gt~=apply_mapping(seg_2_body_map, as_to_gt_map))-1;


% Number of mergers: all segments split from their auto. seg bodies must be
% merged with their correct bodies. We must also merge auto. seg. bodies
% assigned to the same gt. body.

% Number of auto. seg. bodies assigned to each body
[gt_label, count] = count_row_occurence(as_to_gt_map(:,2));
gt_n_as_member_body = [];
gt_n_as_member_body(gt_label+1) = count;

% If N auto. seg. bodies are assigned to one gt. body then we have to make
% N-1 mergers.
n_merge_as_body = sum(nonzeros(gt_n_as_member_body)-1);

% All segments split from their auto. seg. bodies need to be merged to the
% correct body.
n_merge = n_merge_as_body + n_split;

% List of segments that need to be merged into some other body
as_to_gt_map_overlap = sortrows([as_to_gt_overlap, as_to_gt_map], 1);
[junk, I] = unique(as_to_gt_map_overlap(:,3), 'last');
as_bodies_correctly_assigned = as_to_gt_map_overlap(I,2);
merge_segment_id = ...
  nonzeros(find(~ismember(seg_2_body_map, [0; as_bodies_correctly_assigned]))-1);
return
end
