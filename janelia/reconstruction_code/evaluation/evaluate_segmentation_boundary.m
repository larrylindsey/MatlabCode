function [precision, recall] = evaluate_segmentation_boundary(...
  segmentation_map, segmentation_map_gt, eval_parameters, image, ...
  oversegment_ignore_map)
% [precision, recall] = evaluate_segmentation_boundary(segmentation_map,
%   segmentation_map_gt, eval_parameters, image, oversegment_ignore_map)
%
% Evaluate the segmentation boundaries using CSA++ - bipartite matching
% with groundtruth segment boundaries.
% See Learning to Detect Natural Image Boundaries Using Local Brightness,
% Color and Texture Cues - Martin et al. PAMI04.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

boundary_map_gt = double(segmentation_map_gt==0);
n_boundary_points_gt = sum(boundary_map_gt(:));

% apply watershed to ensure that only closed segment boundaries are
% counted as detections
boundary_map = segmentation_map==0;
boundary_map_ws = watershed(boundary_map,4);
boundary_map = double(boundary_map_ws==0);
% figure(880); imshow(boundary_map, []);

% match with groundtruth segmentation map
[match, match_gt] = correspondPixels(boundary_map, boundary_map_gt, ...
  eval_parameters.dmax/sqrt(size(boundary_map,1)^2+size(boundary_map,2)^2));

% remove isolated unmatched pixels - may be due to small variations in
% segment boundaries.
unmatch_gt = match_gt<=0 & boundary_map_gt>0;
unmatch_gt = bwareaopen(unmatch_gt, 4);
match_gt = unmatch_gt<=0 & boundary_map_gt>0;

unmatch = match<=0 & boundary_map>0;
unmatch = bwareaopen(unmatch, 4);
match = unmatch<=0 & boundary_map>0;

% ignore certain over-segmentations if so configured
if((exist('oversegment_ignore_map', 'var')==1) && ~isempty(oversegment_ignore_map))
  match_masked = (match>0) | ((boundary_map==1) & (match==0) & oversegment_ignore_map);
else
  match_masked = match>0;
end

recall = [sum(match_gt(:)>0), n_boundary_points_gt];
precision = [sum(match_masked(:)>0), sum(boundary_map(:))];

if(eval_parameters.save_matching)
  [py, px] = find(boundary_map_gt);
  marker.color = 'g'; marker.type = '.';
  match_map_gt = plot_markers(repmat(image, [1 1 3]), px, py, marker);
  [py, px] = find(boundary_map_gt.*(match_gt==0));
  marker.color = 'm'; marker.type = '.';
  match_map_gt = plot_markers(match_map_gt, px, py, marker);
  figure(881);
  axes('Position', [0 0 1 1]);
  set(gcf, 'Position', [20, 20, size(image, 2), size(image, 1)]);
  imshow(match_map_gt);
  if(~isfield(eval_parameters, 'write_matching') || ...
      eval_parameters.write_matching)
    imwrite(match_map_gt, get_storage_file_name([eval_parameters.save_prefix, ...
      '.gt.tif']));
  end
  
  [py, px] = find(boundary_map);
  marker.color = 'g'; marker.type = '.';
  match_map = plot_markers(repmat(image, [1 1 3]), px, py, marker);
  [py, px] = find(boundary_map.*(match==0));
  marker.color = 'b'; marker.type = '.';
  match_map = plot_markers(match_map, px, py, marker);
  [py, px] = find(boundary_map.*(match_masked==0));
  marker.color = 'r'; marker.type = '.';
  match_map = plot_markers(match_map, px, py, marker);
  figure(882);
  axes('Position', [0 0 1 1]);
  set(gcf, 'Position', [1500, 20, size(image, 2), size(image, 1)]);
  imshow(match_map);
  if(~isfield(eval_parameters, 'write_matching') || ...
      eval_parameters.write_matching)
    imwrite(match_map, get_storage_file_name([eval_parameters.save_prefix, ...
      '.as.tif']));
  end
end;

end
