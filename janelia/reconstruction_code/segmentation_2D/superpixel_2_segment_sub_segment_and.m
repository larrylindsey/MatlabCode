function superpixel_2_segment_sub_segment_and(superpixel_seg, boundary_map, image, ...
  seg_config, version_name, save_prefix, superpixel_suffix, proofread_seg, ...
  config, image_prefix)
% superpixel_2_segment_sub_segment_and(seg_configboundary_map, image, ...
%   seg_config, version_name, save_prefix)
% Perform 2D segmentation based on superpixels by "and"ing two more segment
% maps.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

if(seg_config.is_verbose)
  fprintf('START: superpixel_2_segment_sub_segment_and\n');
end

if(nargin<=7)
  proofread_seg = [];
end

if(isfield(superpixel_seg, 'locked_labels') && ...
    ~isempty(superpixel_seg.locked_labels) && ...
    isempty(proofread_seg))
  error(['The superpixel map has locked labels but I am not being provided ', ...
    'any proofread segment map for correctly mapping these.']);
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

label_map = superpixel_seg.label_map;
superpixel_to_seg_label = [];
superpixel_to_seg_label(:,1) = nonzeros(unique(label_map(:)));
superpixel_to_seg_label(:,2) = superpixel_to_seg_label(:,1);
save_suffix_for_and = [];
fprintf('number of superpixels: %d\n', size(superpixel_to_seg_label,1));
for i = 1:length(seg_config.segment_and.to_merge_segment_maps)
  fprintf('seg_config.segment_and.to_merge_segment_maps(%d).seg_suffix: %s\n', ...
    i, seg_config.segment_and.to_merge_segment_maps(i).seg_suffix);
  sg = load2([get_reconstruction_dir(config), seg_config.dir, ...
    seg_config.segment_and.to_merge_segment_maps(i).method, '/', ...
    image_prefix, ...
    seg_config.segment_and.to_merge_segment_maps(i).seg_suffix, '.mat']);
  [label_map, junk, sp_to_seg] = segment_and(label_map, sg.label_map, ...
    boundary_map<0.5, seg_config.segment_and.overlap_area_threshold); %#ok<ASGLU>
  superpixel_to_seg_label(:,2) = apply_mapping(superpixel_to_seg_label(:,2), ...
    sp_to_seg);
  
  fprintf('number of segments: %d\n', nnz(unique(superpixel_to_seg_label(:,2))));
  
  save_suffix_for_and = [save_suffix_for_and, ...
    seg_config.segment_and.to_merge_segment_maps(i).save_suffix]; %#ok<AGROW>
end

save_file_name_suffix = ['.sa_sp_', num2str(seg_config.segment_and.overlap_area_threshold), ...
  save_suffix_for_and, superpixel_suffix, '.', version_name];

if(seg_config.is_verbose)
  fprintf('save_file_name_suffix: %s\n', save_file_name_suffix);
end

if(isfield(superpixel_seg, 'locked_labels') && ...
    ~isempty(superpixel_seg.locked_labels) && ...
    ~isempty(proofread_seg))
  if(seg_config.is_verbose)
    fprintf('Retaining locked segments ensuring non-overlapping labels\n');
  end
  if(~isempty(label_map))
    % make labels unique
    seg_ids = nonzeros(unique(label_map(:)));
    seg_ids_proofread = nonzeros(unique(proofread_seg.superpixel_to_seg_label(:,2)));
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
  
  superpixel_to_seg_map_proofread = [];
  superpixel_to_seg_map_proofread(1+proofread_seg.superpixel_to_seg_label(:,1)) = ...
    proofread_seg.superpixel_to_seg_label(:,2);
  superpixel_label_map_proofread = superpixel_seg.original_label_map;
  superpixel_label_map_proofread(~superpixel_seg.locked_mask) = 0;
  % some small pixels might have been removed during the locked mask
  % erosion and dilation
  superpixel_label_map_proofread(~ismember(superpixel_label_map_proofread, ...
    proofread_seg.superpixel_to_seg_label(:,1))) = 0;
  segment_label_map_proofread = superpixel_to_seg_map_proofread(1+ ...
    superpixel_label_map_proofread);
  label_map(superpixel_seg.locked_mask) = ...
    segment_label_map_proofread(superpixel_seg.locked_mask);
  
  superpixel_to_seg_label = [superpixel_to_seg_label; ...
    proofread_seg.superpixel_to_seg_label];
  locked_labels = proofread_seg.locked_labels;
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
  fprintf('STOP: superpixel_2_segment_sub_segment_and\n');
end
return
end
