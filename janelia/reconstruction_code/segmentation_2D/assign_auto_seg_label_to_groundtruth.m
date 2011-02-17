function assignment = assign_auto_seg_label_to_groundtruth(segment_map_as, ...
  segment_map_gt, param)
% assignment = assign_auto_seg_label_to_groundtruth(segment_map_as,
%   segment_map_gt, param)
% param           Parameters for including automatic segments that overlap
%                   with multiple ground-truth segments.
% param.normalized_area_threshold   minimum value of overlap area / area of
%                                     automatic segment. 0.9 seems
%                                     reasonable.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%

% Get amount of overlap between automatic and ground-truth segments in
% number of pixels
label_pairs = [segment_map_as(:), segment_map_gt(:)];
label_pairs = label_pairs(label_pairs(:,1)>0,:);
[label_pairs, label_pairs_n_pixel] = count_row_occurence(label_pairs);

% Find automatic segments overlapping with multiple ground-truth segments.
% Automatic segments overlapping with a unique ground-truth segment are
% definitely assigned. Those automatic segments that overlap with more than
% one ground-truth are evaluated for inclusion in the assignment to see if
% assignment is unambiguous.
[labels_as, n_overlap_label_gt] = count_row_occurence(label_pairs(:,1));
labels_as_multi_overlap = labels_as(n_overlap_label_gt>1);

% Decide whether to generate assignments for automatic segments overlapping
% with multiple ground-truth segments. Quantify ambiguity.
[l_id, counts] = count_row_occurence(segment_map_as(:));
area_as = [];
area_as(1+l_id) = counts;
area_as = area_as';

label_pairs_normalized_area = label_pairs_n_pixel./area_as(1+label_pairs(:,1));
label_pairs_norm_area = [label_pairs, label_pairs_normalized_area];
label_pairs_norm_area = label_pairs_norm_area(...
  label_pairs_norm_area(:,3)>param.normalized_area_threshold,:);
label_pairs_norm_area_s = sortrows(label_pairs_norm_area, 3);
label_pairs_norm_area_s = label_pairs_norm_area_s(...
  ismember(label_pairs_norm_area_s(:,1), labels_as_multi_overlap), :);
[junk, I] = unique(label_pairs_norm_area_s(:,1), 'last');
label_pairs_multi_overlap_to_retain = label_pairs_norm_area_s(I,1:2);

assignment = [label_pairs(~ismember(label_pairs(:,1), labels_as_multi_overlap), :); ...
  label_pairs_multi_overlap_to_retain];

return
end
