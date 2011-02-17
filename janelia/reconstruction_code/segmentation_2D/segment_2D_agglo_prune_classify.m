function [label_map_new, sp_to_seg_map, segment_suffix] = ...
  segment_2D_agglo_prune_classify(label_map, boundary_map, merge_criterion, ...
  merge_criterion_param, image, config, image_prefix, is_verbose)
% [label_map_new, sp_to_seg_map] = segment_2D_agglo_prune_classify(...
%   label_map, boundary_map, merge_criterion, merge_criterion_param)
%
% Shiv Vitaladevuni, Shaul Druckmann
% Janelia Farm Research Campus, HHMI.
%

if(nargin<8)
  is_verbose = true;
end

if(is_verbose)
  fprintf('START: segment_2D_agglo_prune_classify\n');
end
label_map = double(label_map);
boundary_map = double(boundary_map);


% get adjacent label pairs and boundary coordinates
boundary_sets = get_segment_boundary_sets(label_map, boundary_map);


% call the merge_criterion
to_merge_sets = [];
label_map_new = label_map;
sp_ids = [0; nonzeros(unique(label_map(:)))];
sp_to_seg_map = [sp_ids, sp_ids];
switch(merge_criterion)
  case 'min'
    if(max(label_map(:))>0)
      to_merge_sets = get_merge_decision_min(...
        label_map, boundary_map, boundary_sets, merge_criterion_param);
    else
      to_merge_sets = [];
    end
    segment_suffix = [merge_criterion, '_', ...
      num2str(merge_criterion_param.min_threshold)];
  case 'region_max'
    if(max(label_map(:))>0)
      to_merge_sets = get_merge_decision_region_max(...
        label_map, boundary_map, boundary_sets, merge_criterion_param);
    else
      to_merge_sets = [];
    end
    segment_suffix = [merge_criterion, '_', ...
      num2str(merge_criterion_param.patch_size), '_', ...
      num2str(merge_criterion_param.min_threshold)];
  case 'boundary_hist_boost'
    if(max(label_map(:))>0)
      [label_map_new, sp_to_seg_map] = get_merge_decision_boundary_hist_boost(...
        label_map, boundary_map, boundary_sets, merge_criterion_param, image, config);
      to_merge_sets = [];
    end
    segment_suffix = [merge_criterion, '_', ...
      num2str(merge_criterion_param.threshold, '%g')];
  case 'boundary_hist_wekarf'
    if(max(label_map(:))>0)
      [label_map_new, sp_to_seg_map] = get_merge_decision_boundary_hist_wekarf(...
        label_map, boundary_map, boundary_sets, merge_criterion_param, image, config);
      to_merge_sets = [];
    end
    segment_suffix = [merge_criterion, '_', ...
      num2str(merge_criterion_param.threshold, '%g')];
  case 'merge_mito_segment'
    if(max(label_map(:))>0)
      to_merge_sets = get_merge_decision_mito_segment(...
        label_map, image, boundary_sets, merge_criterion_param, image_prefix, config);
    else
      to_merge_sets = [];
    end
    segment_suffix = [merge_criterion, '_', ...
      num2str(1), '_', ...
      num2str(0.9)];
  otherwise
    error('prune classify merge criterion not recognized.');
end

if(~isempty(to_merge_sets))
  % make mergers
  max_label = max(label_map(:));
  % merger adjacency matrix for connected components analysis
  A = sparse([], [], [], max_label, max_label, 2*size(to_merge_sets,1));
  A(sub2ind(size(A), to_merge_sets(:,1), to_merge_sets(:,2))) = 1;
  A(sub2ind(size(A), to_merge_sets(:,2), to_merge_sets(:,1))) = 1;
  
  seg_label_mapping = [0; get_connected_components_c(A, 0.5)];
  
  % generate new label map
  label_map_new = seg_label_mapping(1+label_map);
  label_map_new = double(remove_merged_boundaries_2D(uint32(label_map_new)));
  
  % generate superpixel to segment map
  sp_to_seg_map = [find(seg_label_mapping)-1, ...
    seg_label_mapping(seg_label_mapping>0)];
end

if(is_verbose)
  fprintf('STOP: segment_2D_agglo_prune_classify\n');
end
return
end
