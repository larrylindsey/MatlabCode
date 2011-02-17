function segment_2_seg_sub_agglo_AMB_ABIFL_AMdB(...
  superpixel_seg, boundary_map, image, seg_config, version_name, save_prefix, ...
  superpixel_suffix, exclusion_mask)
% superpixel_2_seg_grayscale_agglo_boundary_vs_interior_FisherLD(...
%   superpixel_seg, boundary_map, image, seg_config, version_name, save_prefix, ...
%   superpixel_suffix, exclusion_mask)
% Perform agglomerative (merger only) segmentation starting from
% superpixels. Criterion: Merge two segments if Fisher Linear Discriminant
% between boundary values and segment-interior values is less than a
% threshold.
%
% Part of iterative agglomerative clustering methods
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v0  08122008  init code
%

if(seg_config.is_verbose)
  fprintf('START: segment_2_seg_sub_agglo_AMB_ABIFL_AMdB\n');
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

previous_segment_map = superpixel_seg.label_map;

superpixel_2_seg_map = zeros(max(superpixel_seg.superpixel_to_seg_label(:,1))+1,1);
superpixel_2_seg_map(superpixel_seg.superpixel_to_seg_label(:,1)+1) = ...
  superpixel_seg.superpixel_to_seg_label(:,2);
final_superpixel_2_seg_map = superpixel_2_seg_map;

switch(seg_config.method)
  case 'grayscale_AglBIFL'
    if(seg_config.is_verbose)
      fprintf('computing agglo. boundary vs. interior FL ...\n');
    end
    method_save_suffix = 'gs_abif';
  case 'grayscale_AglMeanB'
    if(seg_config.is_verbose)
      fprintf('computing agglo. mean boundary ...\n');
    end
    method_save_suffix = 'gs_amb';
  case 'grayscale_AglMedianB'
    if(seg_config.is_verbose)
      fprintf('computing agglo. median boundary ...\n');
    end
    method_save_suffix = 'gs_amdb';
  otherwise
    error('Segmentation method not recognized.');
end

f_threshold = 0;
for t = 1:length(seg_config.f_threshold_seq)
  f_threshold_prev = f_threshold;
  f_threshold = seg_config.f_threshold_seq(t);
  length_threshold = seg_config.length_threshold_seq(t);
  save_file_name_suffix = ['.', method_save_suffix, '_sp_T', num2str(f_threshold), ...
    '_', num2str(f_threshold_prev), '_b', num2str(length_threshold), ...
    superpixel_suffix, '.', version_name];
  
  if(seg_config.is_verbose)
    fprintf('f_threshold: %g, length_threshold: %g\n', f_threshold, length_threshold);
    fprintf('save_file_name_suffix:%s\n', save_file_name_suffix);
  end
  
  if(~is_blank_segmentation)
    segment_boundary_d = imdilate(previous_segment_map == 0, strel('disk', 5));
    exclusion_mask_current = max(exclusion_mask, segment_boundary_d);
    
    switch(seg_config.method)
      case 'grayscale_AglBIFL'
        [label_map, superpixel_to_seg_label] = ...
          seg_frm_superpixel_boundary_vs_interior_FisherLD(previous_segment_map, ...
          boundary_map, f_threshold, length_threshold, seg_config.max_area_threshold, ...
          exclusion_mask_current);
      case 'grayscale_AglMeanB'
        [label_map, superpixel_to_seg_label] = ...
          segmentation_from_superpixels_mean_boundary_value(previous_segment_map, ...
          boundary_map, f_threshold, length_threshold);
      case 'grayscale_AglMedianB'
        [label_map, superpixel_to_seg_label] = ...
          segmentation_from_superpixels_median_boundary_value(previous_segment_map, ...
          boundary_map, f_threshold, length_threshold);
    end
    
    superpixel_2_seg_map = zeros(max(superpixel_to_seg_label(:,1))+1,1);
    superpixel_2_seg_map(superpixel_to_seg_label(:,1)+1) = ...
      superpixel_to_seg_label(:,2);
    final_superpixel_2_seg_map = superpixel_2_seg_map(final_superpixel_2_seg_map+1);
    
    superpixel_to_seg_label = [find(final_superpixel_2_seg_map)-1; ...
      final_superpixel_2_seg_map(final_superpixel_2_seg_map>0)]';
  else
    label_map = zeros(size(superpixel_seg.label_map));
    superpixel_to_seg_label = [];
  end
  previous_segment_map = label_map;
  
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
        seg_id_map(1+seg_ids) = seg_ids; %#ok<AGROW>
        seg_id_map(1+seg_ids(non_unique)) = possible_ids(1:n_non_unique); %#ok<AGROW>
        label_map = seg_id_map(1+label_map);
        superpixel_to_seg_label(:,2) = seg_id_map(1+superpixel_to_seg_label(:,2));
      end
    end
    
    label_map(superpixel_seg.locked_mask) = ...
      superpixel_seg.original_label_map(superpixel_seg.locked_mask);
    
    superpixel_to_seg_label = [superpixel_to_seg_label; ...
      superpixel_seg.original_superpixel_to_seg_label(ismember( ...
      superpixel_seg.original_superpixel_to_seg_label(:,2), ...
      superpixel_seg.locked_labels), :)]; %#ok<AGROW>
    
    locked_labels = superpixel_seg.locked_labels;
  else
    locked_labels = [];
  end
  
  label_map(label_map<0) = 0;
  
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

  % retain only those segments that are present in the raster label map
  seg_ids = unique(label_map(:));
  if(~isempty(superpixel_to_seg_label))
    superpixel_to_seg_label = superpixel_to_seg_label(...
      ismember(superpixel_to_seg_label(:,2), seg_ids), :);
  end
  if(~isempty(locked_labels))
    locked_labels = locked_labels(ismember(locked_labels, seg_ids)); %#ok<NASGU>
  end
  
  if(~isfield(seg_config, 'save_f_thresholds') || ...
      ismember(single(f_threshold), single(seg_config.save_f_thresholds))==1)
    save2([save_prefix, save_file_name_suffix,'.mat'], ...
      'label_map', 'superpixel_to_seg_label', 'locked_labels');
    
    if(isfield(seg_config, 'save_segmentation_overlay') && ...
        seg_config.save_segmentation_overlay)
      segment_boundary = label_map==0;
      save_image = [];
      save_image(:,:,1) = image;
      save_image(:,:,2) = image+segment_boundary;
      save_image(:,:,3) = image;
      imwrite(save_image, ...
        get_storage_file_name([save_prefix, save_file_name_suffix,'.tif']));
    end
  end;
end;
if(seg_config.is_verbose)
  fprintf('done.\n');
end
if(seg_config.is_verbose)
  fprintf('STOP: segment_2_seg_sub_agglo_AMB_ABIFL_AMdB\n');
end
return;
end
