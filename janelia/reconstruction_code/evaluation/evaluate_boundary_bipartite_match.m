function match_result = evaluate_boundary_bipartite_match(boundary_map_as, boundary_map_gt, ...
  eval_parameters, oversegment_ignore_map)
% match_result = evaluate_boundary_bipartite_match(boundary_map_as, boundary_map_gt, ...
%   eval_parameters, oversegment_ignore_map)

boundary_map_as = double(boundary_map_as);
boundary_map_gt = double(boundary_map_gt);

n_boundary_points_as = sum(boundary_map_as(:));
n_boundary_points_gt = sum(boundary_map_gt(:));

% match with groundtruth segmentation map
[match_as, match_gt] = correspondPixels(boundary_map_as, boundary_map_gt, ...
  eval_parameters.dmax/sqrt(size(boundary_map_as,1)^2+size(boundary_map_as,2)^2));

% remove isolated unmatched pixels - may be due to small variations in
% segment boundaries.
unmatch_gt = match_gt<=0 & boundary_map_gt>0;
unmatch_gt = bwareaopen(unmatch_gt, 4);
is_matched_gt = unmatch_gt<=0 & boundary_map_gt>0;

unmatch_as = match_as<=0 & boundary_map_as>0;
unmatch_as = bwareaopen(unmatch_as, 4);
is_matched_as = unmatch_as<=0 & boundary_map_as>0;

% ignore certain over-segmentations if so configured
if((exist('oversegment_ignore_map', 'var')==1) && ~isempty(oversegment_ignore_map))
  is_matched_as_masked = (is_matched_as>0) | ((boundary_map_as==1) & (is_matched_as==0) & oversegment_ignore_map);
else
  is_matched_as_masked = is_matched_as;
end

match_result = [];
match_result.recall = [sum(is_matched_gt(:)>0), n_boundary_points_gt];
match_result.precision = [sum(is_matched_as_masked(:)>0), n_boundary_points_as];
match_result.match_as = match_as;
match_result.match_as(match_as==0 & boundary_map_as>0) = -1;
match_result.match_gt = match_gt;
match_result.match_gt(match_gt==0 & boundary_map_gt>0) = -1;
return;
end
