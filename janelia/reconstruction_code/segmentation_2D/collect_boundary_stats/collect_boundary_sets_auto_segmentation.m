function seg_boundary_sets = collect_boundary_sets_auto_segmentation(...
  segmentation_map, segmentation_map_gt, eval_parameters, image_boundary_map, ...
  oversegment_ignore_map)
% seg_boundary_sets = collect_boundary_sets_auto_segmentation(...
%  segmentation_map, segmentation_map_gt, eval_parameters, image_boundary_map, ...
%  oversegment_ignore_map)
% Collect sets of segment boundary elements from automatic segmentations.
% Compare the segmentation with groundtruth and divide into two classes:
% correct and over-segmentations.
%
% For bipartite matching see Learning to Detect Natural Image Boundaries
% Using Local Brightness, Color and Texture Cues - Martin et al. PAMI04.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

boundary_map_gt = double(segmentation_map_gt==0);

% % apply watershed to ensure that only closed segment boundaries are
% % counted as detections
% boundary_map = segmentation_map==0;
% boundary_map_ws = watershed(boundary_map,4);
% boundary_map = double(boundary_map_ws==0);
boundary_map = double(segmentation_map==0);
figure(880); imshow(boundary_map, []);

% match with groundtruth segmentation map
[match, match_gt] = correspondPixels(boundary_map, boundary_map_gt, ...
  eval_parameters.dmax/sqrt(size(boundary_map,1)^2+size(boundary_map,2)^2)); %#ok<NASGU>

% ignore certain over-segmentations if so configured
if((exist('oversegment_ignore_map', 'var')==1) && ~isempty(oversegment_ignore_map))
  seg_boundary_matched = (match>0) & (oversegment_ignore_map==0);
  seg_boundary_unmatched = (boundary_map>0) & (match==0) & (oversegment_ignore_map==0);
else
  seg_boundary_matched = match>0;
  seg_boundary_unmatched = (boundary_map>0) & (match==0);
end

% remove small clusters of (2/3) boundary elements that are labelled as
% unmatched from the unmatched mask and put in matched mask. These are
% erroneous error flags cause due to one-to-one matching
unmatched_significant = bwareaopen(seg_boundary_unmatched, 3);
seg_boundary_matched_add = (1-unmatched_significant).*seg_boundary_unmatched;
seg_boundary_matched = max(seg_boundary_matched, seg_boundary_matched_add);

seg_boundary_matched = double(seg_boundary_matched);
seg_boundary_unmatched = double(seg_boundary_unmatched);
% figure(1); imshow(boundary_map_gt);
% figure(2); imshow(boundary_map);
% figure(3); imshow(seg_boundary_matched);
% figure(4); imshow(seg_boundary_unmatched);

seg_boundary_sets.matched = get_segment_boundary_sets_masked(segmentation_map, ...
  image_boundary_map, seg_boundary_matched);
seg_boundary_sets.unmatched = get_segment_boundary_sets_masked(segmentation_map, ...
  image_boundary_map, seg_boundary_unmatched);

return
end
