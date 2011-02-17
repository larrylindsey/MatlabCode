function superpixel_2_superpixel_sub_prune_classify(superpixel_seg, boundary_map, image, ...
  seg_config, version_name, save_prefix, superpixel_suffix)
% superpixel_2_superpixel_sub_prune_classify(seg_configboundary_map, image, ...
%   seg_config, version_name, save_prefix)
% Perform 2D segmentation based on superpixels using agglomerative
% prune-classify
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

if(seg_config.is_verbose)
  fprintf('START: superpixel_2_superpixel_sub_prune_classify\n');
end

superpixel_seg.original_label_map = superpixel_seg.label_map;
if(isfield(superpixel_seg, 'locked_labels') && ...
    ~isempty(superpixel_seg.locked_labels))
  superpixel_seg.locked_mask = ismember(superpixel_seg.label_map, ...
    superpixel_seg.locked_labels);
  superpixel_seg.locked_mask = imerode(imdilate(superpixel_seg.locked_mask, ...
    strel('disk', 4)), strel('disk', 4));
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

[label_map, junk, seg_config.segment_suffix] = ...
  segment_2D_agglo_prune_classify(superpixel_seg.label_map, ...
  boundary_map, seg_config.prune_classify.merge_criterion, ...
  seg_config.prune_classify.merge_criterion_param); %#ok<ASGLU>

save_file_name_suffix = ['.apc_sp_', seg_config.segment_suffix, ...
  superpixel_suffix, '.', version_name];

if(seg_config.is_verbose)
  fprintf('save_file_name_suffix: %s\n', save_file_name_suffix);
  fprintf('segment_suffix: .apc_sp_%s\n', seg_config.segment_suffix);
end

if(isfield(superpixel_seg, 'locked_labels') && ...
    ~isempty(superpixel_seg.locked_labels))
  if(seg_config.is_verbose_figures)
    fprintf('Retaining locked superpixels ensuring non-overlapping labels\n');
  end
  if(~isempty(label_map))
    % make labels unique
    sp_ids = nonzeros(unique(label_map(:)));
    non_unique = ismember(sp_ids, superpixel_seg.locked_labels);
    n_non_unique = nnz(non_unique);
    if(n_non_unique>0)
      possible_ids = 1:(max([superpixel_seg.locked_labels; sp_ids])+n_non_unique);
      possible_ids = possible_ids(~ismember(possible_ids, ...
        [superpixel_seg.locked_labels; sp_ids]));
      sp_id_map = [];
      sp_id_map(1+sp_ids) = sp_ids;
      sp_id_map(1+sp_ids(non_unique)) = possible_ids(1:n_non_unique);
      label_map = sp_id_map(1+label_map);
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

save2([save_prefix, save_file_name_suffix,'.mat'], ...
  'label_map', 'locked_labels');

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

if(seg_config.is_verbose)
  fprintf('STOP: superpixel_2_superpixel_sub_prune_classify\n');
end
return
end
