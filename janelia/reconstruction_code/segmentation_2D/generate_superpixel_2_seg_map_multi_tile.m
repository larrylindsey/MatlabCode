function generate_superpixel_2_seg_map_multi_tile(config)
% generate_superpixel_2_seg_map_multi_tile(config)
% generate superpixel_2_seg_map from superpixel_to_seg_label obtained from
% segmentation. Multiple tiles per section (trakEM xml)
%
% superpixel_2_seg_map{z} is a Nx1 matrix, where N is the number of
% superpixels in the image.
% superpixel_2_seg_map{z}(s) = id. of segment containing the superpixel.
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
% v1  09302008  modified for multiple tiles per section (trakEM xml)
%

fprintf('Constructing superpixel_2_seg_map ..\n');
stack_config = config.stack;
seg_config = config.superpixel_2_seg(end);
linkage_config = config.linkage;
apply_config = linkage_config.apply;

seg_dir = [get_reconstruction_dir(config), seg_config.dir, seg_config.method, '/'];

superpixel_2_seg_map = [];

case_id = stack_config.case_ids(1);
fprintf('%d ', case_id);

segment_2D_label_map_1 = load(sprintf([seg_dir, stack_config.image_prefix, ...
  apply_config.segmentation_suffix,'.mat'] , case_id));
n_superpixel = max(segment_2D_label_map_1.superpixel_to_seg_label(:,1));
superpixel_2_seg_map{end+1} = zeros(n_superpixel+1, 1);
superpixel_2_seg_map{end}(segment_2D_label_map_1.superpixel_to_seg_label(:,1)+1) = ...
  segment_2D_label_map_1.superpixel_to_seg_label(:,2);

for i = 2:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  fprintf('%d ', case_id);
  if(mod(i,20)==0)
    fprintf('\n');
  end;
  
  segment_2D_label_map_2 = load(sprintf([seg_dir, stack_config.image_prefix, ...
    apply_config.segmentation_suffix,'.mat'] , case_id));
  n_superpixel = max(segment_2D_label_map_2.superpixel_to_seg_label(:,1));
  superpixel_2_seg_map{end+1} = zeros(n_superpixel+1, 1);
  superpixel_2_seg_map{end}(segment_2D_label_map_2.superpixel_to_seg_label(:,1)+1) = ...
    segment_2D_label_map_2.superpixel_to_seg_label(:,2);
end;


save([seg_dir, 'superpixel_2_seg_map', apply_config.segmentation_suffix, '.mat'], 'superpixel_2_seg_map');

fprintf('\ndone.\n');

