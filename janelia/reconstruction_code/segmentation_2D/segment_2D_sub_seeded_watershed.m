function segment_2D_sub_seeded_watershed(segment_map_init, boundary_map, image, ...
  fold_mask, watershed_seed_mask, seg_config, version_name, save_prefix, proofread_sp)
% segment_2D_sub_seeded_watershed(segment_map_init, boundary_map, image, ...
%   seg_config, version_name, save_prefix)
% Perform 2D segmentation with seeded watershed
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v0  ~Feb.2008 init code
% v1  04112008  modified for reconstruction pipeline
% v2  02112009  modified to unify common code among algorithms
%

if(seg_config.is_verbose)
  fprintf('START: segment_2D_sub_seeded_watershed\n');
end

if(nargin<=6)
  proofread_sp = [];
end

is_blank_segmentation = max(segment_map_init(:))<=1;
if(is_blank_segmentation)
  boundary_map = ones(size(segment_map_init));
end

save_file_name_suffix = ['.sws.', version_name];
if(seg_config.is_verbose)
  fprintf('save_file_name_suffix:%s\n', save_file_name_suffix);
end

if(~is_blank_segmentation)
  label_map = compute_seeded_watershed_segment_and(...
    segment_map_init, boundary_map, watershed_seed_mask, fold_mask, seg_config);
  label_map = double(remove_merged_boundaries_2D(uint32(label_map)));
else
  label_map = zeros(size(segment_map_init));
end

if(~isempty(proofread_sp))
  if(seg_config.is_verbose)
    fprintf('Combining new seg. with proofread seg. ensuring non-overlapping labels\n');
  end
  if(~isempty(label_map))
    % make labels unique
    sp_ids = unique(label_map(:));
    non_unique = ismember(sp_ids, proofread_sp.locked_labels);
    n_non_unique = nnz(non_unique);
    if(n_non_unique>0)
      possible_ids = 1:(max([proofread_sp.locked_labels; sp_ids])+n_non_unique);
      possible_ids = possible_ids(~ismember(possible_ids, ...
        [proofread_sp.locked_labels; sp_ids]));
      sp_id_map = [];
      sp_id_map(1+sp_ids) = sp_ids;
      sp_id_map(1+sp_ids(non_unique)) = possible_ids(1:n_non_unique);
      label_map = sp_id_map(1+label_map);
    end
  end
  label_map(proofread_sp.locked_mask) = ...
    proofread_sp.label_map(proofread_sp.locked_mask);
  locked_labels = proofread_sp.locked_labels; %#ok<NASGU>
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

if(~isfield(seg_config, 'save_f_thresholds') || ...
    ismember(single(f_threshold), single(seg_config.save_f_thresholds))==1)
  save2([save_prefix, save_file_name_suffix,'.mat'],  'label_map', ...
    'locked_labels');
  
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
end

if(seg_config.is_verbose)
  fprintf('STOP: segment_2D_sub_seeded_watershed\n');
end
return;
end
