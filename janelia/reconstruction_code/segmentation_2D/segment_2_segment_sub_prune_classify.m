function segment_2_segment_sub_prune_classify(superpixel_seg, boundary_map, image, ...
  seg_config, version_name, save_prefix, superpixel_suffix, config, image_prefix)
% segment_2_segment_sub_prune_classify(seg_configboundary_map, image, ...
%   seg_config, version_name, save_prefix)
% Perform 2D segmentation based on superpixels using agglomerative
% prune-classify
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

if(seg_config.is_verbose)
  fprintf('START: segment_2_segment_sub_prune_classify\n');
end

superpixel_seg.original_label_map = superpixel_seg.label_map;
superpixel_seg.original_superpixel_to_seg_label = ...
  superpixel_seg.superpixel_to_seg_label;
if(isfield(superpixel_seg, 'locked_labels') && ...
    ~isempty(superpixel_seg.locked_labels))
  superpixel_seg.locked_mask = ismember(superpixel_seg.label_map, ...
    superpixel_seg.locked_labels);
  superpixel_seg.label_map(ismember(superpixel_seg.label_map, ...
    superpixel_seg.locked_labels)) = 0;
  superpixel_seg.superpixel_to_seg_label = ...
    superpixel_seg.superpixel_to_seg_label(~ismember( ...
    superpixel_seg.superpixel_to_seg_label(:,2), ...
    superpixel_seg.locked_labels), :);
end

is_blank_segmentation = max(superpixel_seg.label_map(:))<=1;
if(is_blank_segmentation)
  boundary_map = ones(size(superpixel_seg.label_map));
end

if(is_blank_segmentation)
  superpixel_seg.label_map = zeros(size(superpixel_seg.label_map));
end

[label_map, superpixel_to_seg_label, seg_config.segment_suffix] = ...
  segment_2D_agglo_prune_classify(...
  superpixel_seg.label_map, boundary_map, ...
  seg_config.prune_classify.merge_criterion, ...
  seg_config.prune_classify.merge_criterion_param, image, config, image_prefix);

save_file_name_suffix = ['.apc_sp_', seg_config.segment_suffix, ...
  superpixel_suffix, '.', version_name];

if(seg_config.is_verbose)
  fprintf('segment_suffix: %s\n', seg_config.segment_suffix);
  fprintf('save_file_name_suffix: %s\n', save_file_name_suffix);
end

if(~is_blank_segmentation)
  % carry over the previous stage's superpixel-to-segment mapping
  superpixel_to_seg_map_previous = [];
  superpixel_to_seg_map_previous(1+superpixel_seg.superpixel_to_seg_label(:,1)) = ...
    superpixel_seg.superpixel_to_seg_label(:,2);
  
  superpixel_to_seg_map = [];
  superpixel_to_seg_map(1+superpixel_to_seg_label(:,1)) = ...
    superpixel_to_seg_label(:,2);
  
  superpixel_to_seg_map = superpixel_to_seg_map(1+superpixel_to_seg_map_previous);
  superpixel_to_seg_label = [find(superpixel_to_seg_map)-1; ...
    superpixel_to_seg_map(superpixel_to_seg_map>0)]';
else
  label_map = zeros(size(superpixel_seg.label_map));
  superpixel_to_seg_label = [];
end

if(isfield(superpixel_seg, 'locked_labels') && ...
    ~isempty(superpixel_seg.locked_labels))
  if(seg_config.is_verbose)
    fprintf('Retaining locked segments ensuring non-overlapping labels\n');
  end
  if(~isempty(label_map))
    % make labels unique
    seg_ids = unique(label_map(:));
    seg_ids_proofread = unique(superpixel_seg.locked_labels);
    non_unique = ismember(seg_ids, seg_ids_proofread);
    n_non_unique = nnz(non_unique);
    if(n_non_unique>0)
      possible_ids = 1:(max([seg_ids_proofread; seg_ids])+n_non_unique);
      possible_ids = possible_ids(~ismember(possible_ids, ...
        [seg_ids_proofread; seg_ids]));
      seg_id_map = [];
      seg_id_map(1+seg_ids) = seg_ids;
      seg_id_map(1+seg_ids(non_unique)) = possible_ids(1:n_non_unique);
      label_map = seg_id_map(1+label_map);
      superpixel_to_seg_label(:,2) = seg_id_map(1+superpixel_to_seg_label(:,2));
    end
  end
  
  label_map(superpixel_seg.locked_mask) = ...
    superpixel_seg.original_label_map(superpixel_seg.locked_mask);
  
  superpixel_to_seg_label = [superpixel_to_seg_label; ...
    superpixel_seg.original_superpixel_to_seg_label(ismember( ...
    superpixel_seg.original_superpixel_to_seg_label(:,2), ...
    superpixel_seg.locked_labels), :)];
  
  locked_labels = superpixel_seg.locked_labels;
else
  locked_labels = [];
end

if(seg_config.is_verbose_figures)
  segment_boundary = double(label_map==0);
  [py, px] = find(segment_boundary);
  figure(3);
  imshow(image);
  hold on;
  plot(px, py, '.');
  hold off;
  title('Segmentation boundaries');
end

label_map(label_map<0) = 0;

% retain only those segments that are present in the raster label map
seg_ids = unique(label_map(:));
if(~isempty(superpixel_to_seg_label))
  superpixel_to_seg_label = superpixel_to_seg_label(...
    ismember(superpixel_to_seg_label(:,2), seg_ids), :); %#ok<NASGU>
end
if(~isempty(locked_labels))
  locked_labels = locked_labels(ismember(locked_labels, seg_ids)); %#ok<NASGU>
end

save2([save_prefix, save_file_name_suffix,'.mat'], ...
  'label_map', 'superpixel_to_seg_label', 'locked_labels');
  
if(isfield(seg_config, 'save_segmentation_overlay') && ...
    seg_config.save_segmentation_overlay)
  segment_boundary = label_map==0;
  save_image = [];
  save_image(:,:,1) = image;
  save_image(:,:,2) = image+segment_boundary;
  save_image(:,:,3) = image;
  fn = [save_prefix, save_file_name_suffix,'.tif'];
  imwrite(save_image, get_storage_file_name(fn), 'Description', fn);
end

if(seg_config.is_verbose)
  fprintf('STOP: segment_2_segment_sub_prune_classify\n');
end
return
end
