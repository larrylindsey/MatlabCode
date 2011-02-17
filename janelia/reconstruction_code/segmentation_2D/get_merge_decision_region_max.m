function to_merge_sets = get_merge_decision_region_max(label_map, boundary_map, ...
  boundary_sets, merge_criterion_param) %#ok<INUSL>
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

to_merge_sets = [];
PatchSize = merge_criterion_param.patch_size;
s = size(boundary_map);
for i = 1:length(boundary_sets)
  BorderInd = sub2ind(s,boundary_sets(i).boundary_set(2,:),boundary_sets(i).boundary_set(1,:));
  PatchInd = GetPatchInd(BorderInd',PatchSize,s(1));
  PatchMat = ones(size(PatchInd));
  T = (PatchInd>0 & PatchInd<(s(1).^2));
  PatchMat(T) = boundary_map(PatchInd(T));
  % Workaround for indices outside barrier values
  if (min(max(PatchMat)) < merge_criterion_param.min_threshold)
    to_merge_sets(end+1,:) = ...
      [boundary_sets(i).segment_label_0, boundary_sets(i).segment_label_1]; %#ok<AGROW>
  end
end

return
end
