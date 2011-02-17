function [label_map_1_new, label_map_2_new, seg_label_mapping_1, seg_label_mapping_2] = ...
  segment_and(label_map_1, label_map_2, mask, overlap_threshold)
% [label_map_1_new, label_map_2_new, seg_label_mapping_1, seg_label_mapping_2] = ...
%   segment_and(label_map_1, label_map_2, mask, overlap_threshold)
%
% Merge two segmentations by the area of overlap between the segments
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%

label_map_2(label_map_2>0) = label_map_2(label_map_2>0) + max(label_map_1(:));

lm1 = label_map_1;
lm1(mask==0)=0;
lm2 = label_map_2;
lm2(mask==0)=0;

lp = [lm1(:), lm2(:)];
lp = lp(lp(:,1)>0 & lp(:,2)>0, :);
[label_pairs, count] = count_row_occurence(lp);

A = sparse([], [], [], max(label_map_2(:)), max(label_map_2(:)), ...
  2*size(label_pairs,1));
A(sub2ind(size(A), label_pairs(:,1), label_pairs(:,2))) = count;
A(sub2ind(size(A), label_pairs(:,2), label_pairs(:,1))) = count;

seg_label_mapping = [0; get_connected_components(A, overlap_threshold)];

label_map_1_new = seg_label_mapping(1+label_map_1);
label_map_1_new = double(remove_merged_boundaries_2D(uint32(label_map_1_new)));
seg_label_mapping_1(:,1) = nonzeros(unique(label_map_1(:)));
seg_label_mapping_1(:,2) = seg_label_mapping(1+seg_label_mapping_1(:,1));

label_map_2_new = seg_label_mapping(1+label_map_2);
label_map_2_new = double(remove_merged_boundaries_2D(uint32(label_map_2_new)));
seg_label_mapping_2(:,1) = nonzeros(unique(label_map_2(:)));
seg_label_mapping_2(:,2) = seg_label_mapping(1+seg_label_mapping_2(:,1));
seg_label_mapping_2(:,1) = seg_label_mapping_2(:,1) - max(label_map_1(:));
return
end