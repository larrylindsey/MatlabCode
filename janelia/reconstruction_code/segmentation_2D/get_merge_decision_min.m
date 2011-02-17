function to_merge_sets = get_merge_decision_min(label_map, boundary_map, ...
  boundary_sets, merge_criterion_param) %#ok<INUSL>

to_merge_sets = [];
for i = 1:length(boundary_sets)
  if(min(boundary_sets(i).boundary_set(3,:)<merge_criterion_param.min_threshold))
    to_merge_sets(end+1,:) = ...
      [boundary_sets(i).segment_label_0, boundary_sets(i).segment_label_1]; %#ok<AGROW>
  end
end

return
end
