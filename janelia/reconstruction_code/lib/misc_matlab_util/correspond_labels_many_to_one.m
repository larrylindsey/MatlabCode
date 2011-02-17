function mapping = correspond_labels_many_to_one(labels_1, labels_2, ...
  overlap, area_1, area_2, overlap_threshold, overlap_norm_threshold)
% mapping = correspond_labels_many_to_one(labels_1, labels_2, ...
%   overlap, overlap_threshold, overlap_norm_threshold)
% Corresponds labels from labels_1 to labels_2 in many-to-one mapping.
% Input:
%   labels_1 and labels_2: overlapping labels. labels_1(i) overlaps with
%     labels_2(i).
%   overlap: number of overlapping elements. overlap(i) is the number of
%     overlapping elements between labels_1(i) and labels_2(i).
%   area_1 and area_2: areas of labels.
%   overlap_threshold: overlap(i)>overlap_threshold indicates significant
%     overlap between labels_1(i) and labels_2(i).
%   overlap_norm_threshold: overlap(i)/min(area_1(labels_1(i)), ...
%   /area_2(labels_2(i)))>overlap_norm_threshold indicates significant 
%     overlap between labels_1(i) and labels_2(i)
%
% Output:
%   mapping from labels in labels_1 to labels in labels_2.
%     If a label in labels_1 overlaps significantly with more than one
%     label in labels_2 then it is mapped to 0 (null mapping).
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

overlap_norm = overlap./min(area_1(labels_1+1), area_2(labels_2+1));

significant_overlap = overlap > overlap_threshold & ...
  overlap_norm > overlap_norm_threshold;

label_pairs = [labels_1(significant_overlap), labels_2(significant_overlap)];

[label_1_ids, first_pos] = unique(label_pairs(:,1), 'first');
[label_1_ids, last_pos] = unique(label_pairs(:,1), 'last');

is_mapped_pair_id = first_pos(first_pos==last_pos);

mapping = label_pairs(is_mapped_pair_id, :);

return
end
