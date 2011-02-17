function segment_2_segment_sub_shrink_levelset(superpixel_seg, boundary_map, image, ...
  seg_config, version_name, save_prefix, superpixel_suffix)
% segment_2_segment_sub_shrink_levelset(seg_configboundary_map, image, ...
%   seg_config, version_name, save_prefix)
% Shrink 2D segment regions to high confidence interior regions using
% Level sets.
%
% The idea of modifying segmentations using level sets was developed by
% Pratim Ghosh during his internship at Janelia during summer 2009.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

if(seg_config.is_verbose)
  fprintf('START: segment_2_segment_sub_shrink_levelset\n');
end

superpixel_seg.original_label_map = superpixel_seg.label_map;
if(isfield(superpixel_seg, 'locked_labels') && ...
    ~isempty(superpixel_seg.locked_labels))
  superpixel_seg.locked_mask = ismember(superpixel_seg.label_map, ...
    superpixel_seg.locked_labels);
  superpixel_seg.label_map(ismember(superpixel_seg.label_map, ...
    superpixel_seg.locked_labels)) = 0;
end

is_blank_segmentation = max(superpixel_seg.label_map(:))<=1;
if(is_blank_segmentation)
  boundary_map = ones(size(superpixel_seg.label_map));
end

if(is_blank_segmentation)
  superpixel_seg.label_map = zeros(size(superpixel_seg.label_map));
end

label_map = shrink_segment_levelset_ChanVese(...
  superpixel_seg.label_map, 5*255*boundary_map, ...
  0.5, 0.5, 0.5*255*255, 25, 0.0, 255.0, 0.1);
seg_config.segment_suffix = 'mul5_i25';

save_file_name_suffix = ['.skl_sp_', seg_config.segment_suffix, ...
  superpixel_suffix, '.', version_name];

if(seg_config.is_verbose)
  fprintf('segment_suffix: %s\n', seg_config.segment_suffix);
  fprintf('save_file_name_suffix: %s\n', save_file_name_suffix);
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
    end
  end
  
  label_map(superpixel_seg.locked_mask) = ...
    superpixel_seg.original_label_map(superpixel_seg.locked_mask);
  
  locked_labels = superpixel_seg.locked_labels; %#ok<NASGU>
else
  locked_labels = []; %#ok<NASGU>
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

superpixel_to_seg_label = superpixel_seg.superpixel_to_seg_label; %#ok<NASGU>

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
  fprintf('STOP: segment_2_segment_sub_shrink_levelset\n');
end
return
end
