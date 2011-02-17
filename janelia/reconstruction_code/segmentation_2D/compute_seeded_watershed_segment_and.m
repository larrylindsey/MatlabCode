function segment_map = compute_seeded_watershed_segment_and(...
  segment_map_init, boundary_map, watershed_seed_mask, fold_mask, seg_config)
% segment_map = compute_seeded_watershed_segment_and(...
%   segment_map_init, boundary_map, watershed_seed_mask, seg_config)
% Perform 2D segmentation with seeded watershed
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%

if(seg_config.is_verbose)
  fprintf('START: compute_seeded_watershed_segment_and\n');
end

seeds = {};
seeds{end+1} = watershed_seed_mask;
seeds{end}(1:end-1, 1:end-1) = seeds{end}(2:end, 2:end);
seeds{end+1} = watershed_seed_mask;
seeds{end}(2:end, 2:end) = seeds{end}(1:end-1, 1:end-1);
seeds{end+1} = watershed_seed_mask;
seeds{end}(2:end, 1:end-1) = seeds{end}(1:end-1, 2:end);
seeds{end+1} = watershed_seed_mask;
seeds{end}(1:end-1, 2:end) = seeds{end}(2:end, 1:end-1);
seeds{end+1} = watershed_seed_mask;
seeds{end}(1:end-1, 2:end) = seeds{end}(1:end-1, 1:end-1);
seeds{end+1} = watershed_seed_mask;
seeds{end}(1:end-1, 2:end) = seeds{end}(2:end, 2:end);

segment_map = segment_map_init;
fprintf('num labels: %d\n', nnz(unique(segment_map(:))));
for i = 1:length(seeds)
  bs = imimposemin(boundary_map, seeds{i}, 4);
  ws_bs = watershed(bs, 4);
  if(~isempty(fold_mask))
    ws_bs(fold_mask==0) = 0;
  end
  
  label_map = compute_segmentation_hierarchy_from_watershed_with_min_area_c(...
    ws_bs, bs, 0.1, 200);
  label_map = double(remove_merged_boundaries_2D(uint32(label_map)));
  
  segment_map = segment_and(segment_map, label_map, boundary_map<0.5, 150);
  fprintf('num labels: %d\n', nnz(unique(segment_map(:))));

  label_map = compute_segmentation_hierarchy_from_watershed_with_min_area_c(...
    ws_bs, bs, 0.0125, 200);
  label_map = double(remove_merged_boundaries_2D(uint32(label_map)));
  
  segment_map = segment_and(segment_map, label_map, boundary_map<0.5, 150);
  fprintf('num labels: %d\n', nnz(unique(segment_map(:))));

  label_map = compute_segmentation_hierarchy_from_watershed_with_min_area_c(...
    ws_bs, bs, 0.1, 225);
  label_map = double(remove_merged_boundaries_2D(uint32(label_map)));
  
  segment_map = segment_and(segment_map, label_map, boundary_map<0.5, 150);
  fprintf('num labels: %d\n', nnz(unique(segment_map(:))));

  label_map = compute_segmentation_hierarchy_from_watershed_with_min_area_c(...
    ws_bs, bs, 0.0125, 225);
  label_map = double(remove_merged_boundaries_2D(uint32(label_map)));
  
  segment_map = segment_and(segment_map, label_map, boundary_map<0.5, 150);
  fprintf('num labels: %d\n', nnz(unique(segment_map(:))));
end

if(seg_config.is_verbose)
  fprintf('STOP: compute_seeded_watershed_segment_and\n');
end
return;
end
