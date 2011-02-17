function generate_superpixel_2_seg_map_uni_tile(config)
% generate_superpixel_2_seg_map_uni_tile(config)
% generate superpixel_2_seg_map from superpixel_to_seg_label obtained from
% segmentation. Single tile per section (non trakEM xml)
%
% superpixel_2_seg_map{z} is a Nx1 matrix, where N is the number of
% superpixels in the image.
% superpixel_2_seg_map{z}(s) = id. of segment containing the superpixel.
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
%

fprintf('Constructing superpixel_2_seg_map ..\n');
stack_config = config.stack;
seg_config = config.superpixel_2_seg(end);

case_id = stack_config.case_ids(1);
fprintf('%d ', case_id);

image_prefixes = get_image_prefixes_subdirs(config, case_id);
image_prefix = image_prefixes{1};
[seg_method, seg_suffix] = get_segmentation_suffix(config, image_prefix);
seg_dir = [get_reconstruction_dir(config), seg_config.dir, seg_method, '/'];
segment_2D_label_map_1 = load2([seg_dir, image_prefix, ...
  seg_suffix,'.mat']);
n_superpixel = max(segment_2D_label_map_1.superpixel_to_seg_label(:,1));
superpixel_2_seg_map = zeros(n_superpixel+1, 1);
superpixel_2_seg_map(segment_2D_label_map_1.superpixel_to_seg_label(:,1)+1) = ...
  segment_2D_label_map_1.superpixel_to_seg_label(:,2);
superpixel_ids_sets = segment_2D_label_map_1.superpixel_to_seg_label(:,1)'; %#ok<NASGU>
superpixel_2_seg_map_t = superpixel_2_seg_map; %#ok<NASGU>
save2([get_region_dir(config), 'superpixel_2_seg_map_t.', num2str(case_id), ...
  config.segmentation_choose.choice.seg_suffix, '.mat'], ...
  'superpixel_2_seg_map_t', 'superpixel_ids_sets');


for i = 2:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  fprintf('%d ', case_id);
  if(mod(i,20)==0)
    fprintf('\n');
  end;
  
  image_prefixes = get_image_prefixes_subdirs(config, case_id);
  image_prefix = image_prefixes{1};
  [seg_method, seg_suffix] = get_segmentation_suffix(config, image_prefix);
  seg_dir = [get_reconstruction_dir(config), seg_config.dir, seg_method, '/'];
  segment_2D_label_map_2 = load2([seg_dir, image_prefix, seg_suffix, '.mat']);
  n_superpixel = max(segment_2D_label_map_2.superpixel_to_seg_label(:,1));
  superpixel_2_seg_map = zeros(n_superpixel+1, 1);
  superpixel_2_seg_map(segment_2D_label_map_2.superpixel_to_seg_label(:,1)+1) = ...
    segment_2D_label_map_2.superpixel_to_seg_label(:,2);
  superpixel_ids_sets = segment_2D_label_map_2.superpixel_to_seg_label(:,1)'; %#ok<NASGU>
  superpixel_2_seg_map_t = superpixel_2_seg_map; %#ok<NASGU>
  save2([get_region_dir(config), 'superpixel_2_seg_map_t.',  num2str(case_id), ...
    config.segmentation_choose.choice.seg_suffix, '.mat'], ...
  'superpixel_2_seg_map_t', 'superpixel_ids_sets');
end;


fprintf('\ndone.\n');
return
end
