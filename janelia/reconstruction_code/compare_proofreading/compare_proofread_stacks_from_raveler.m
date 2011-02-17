function compare_proofread_stacks_from_raveler(...
  section_image_name_format, planes, proofread_dir_a, proofread_dir_b, ...
  diff_stack_file_name, rendering_scale, render_original_image)
% compare_proofread_stacks_from_raveler(...
%   section_image_name_format, planes, proofread_dir_a, proofread_dir_b, ...
%   diff_stack_file_name, rendering_scale, render_original_image)
% Compare two proofreadings of a stack obtained from Raveler.
% The two proofreading sessions must first be exported from Raveler.
%
% Input:
%   section_image_name_format
%     specifies file name format with full path in printf format. E.g., if
%     the images were stored as /tmp/a.01.tif, /tmp/a.02.tif, ...
%     /tmp/a.42.tif, then the section_image_name_format should be
%     /tmp/a.%02d.tif and planes should be 1:42.
%   planes
%     list of the section ids.
%   proofread_dir_a/b
%     full path to directories having the exported data. The program
%     expects the superpixel maps to be stored in
%     proofread_dir_a/sp_maps/, and files superpixel-to-segment.txt
%   diff_stack_file_name
%     rendered diff-stack save file name.
%   rendering_scale
%     the amount by which to scale down the sections for rendering. Can be
%     useful for reducing storage space. Default 1.
%   render_original_image
%     whether to render the original image alongside the diff-stack for
%     easy viewing. Default false.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  06192009  init. code
%

if(nargin<=4)
  rendering_scale = 1;
end
if(nargin<=5)
  render_original_image = false;
end

fprintf('Loading proofread stack a ...\n');
seg_0 = uint32(get_body_label_stack_from_raveler_proofread_data(...
  proofread_dir_a, planes, ...
  'superpixel_2_segment_map_file', 'superpixel-to-segment-map.txt', ...
  'segment_2_body_map_file', 'segment-to-body-map.txt', ...
  'superpixel_map_dir', 'sp_maps/', ...
  'superpixel_map_prefix', 'sp_map.%d'));
fprintf('done.\n');
fprintf('Loading proofread stack b ...\n');
seg_1 = uint32(get_body_label_stack_from_raveler_proofread_data(...
  proofread_dir_b, planes, ...
  'superpixel_2_segment_map_file', 'superpixel-to-segment-map.txt', ...
  'segment_2_body_map_file', 'segment-to-body-map.txt', ...
  'superpixel_map_dir', 'sp_maps/', ...
  'superpixel_map_prefix', 'sp_map.%d'));
fprintf('done.\n');

fprintf('Computing difference between proofread stacks ...\n');
min_overlap_norm_threshold = 0.5;
[matched_pairs, diff_stack]  = ...
  compare_label_stacks_bp2(seg_0, seg_1, min_overlap_norm_threshold);
fprintf('done.\n');

fprintf('Rendering difference stack for viewing ...\n');
if(exist(diff_stack_file_name, 'file')==2)
  delete(diff_stack_file_name);
end
if(render_original_image)
  x_disp = 50;
else
  x_disp = 50;
end

for i = 1:length(planes)
  plane = planes(i);
  fprintf('plane: %d\n', plane);
  a = imread(sprintf(section_image_name_format, plane));
  
  if(render_original_image)
    [height, width] = size(a);
    height_2 = round(height/rendering_scale);
    width_2 = round(width/rendering_scale);
    a = imresize(a, [height_2, width_2]);
    
    image = uint8(zeros([height_2, 2*width_2+x_disp, 3]));
    image(:,1:width_2,:) = repmat(a, [1 1 3]);
    image(:,width_2+x_disp+1:end, :) = repmat(a, [1 1 3]);
    
    d = imresize(diff_stack(:,:,i), [height_2, width_2], 'nearest');
    image(:,width_2+x_disp+1:end, 1) = ...
      max(image(:,width_2+x_disp+1:end, 1), uint8(255*(bitand(d,uint8(1))>0)));
    image(:,width_2+x_disp+1:end, 2) = ...
      max(image(:,width_2+x_disp+1:end, 2), uint8(255*(bitand(d,uint8(2))>0)));
    image(:,width_2+x_disp+1:end, 3) = ...
      max(image(:,width_2+x_disp+1:end, 3), uint8(255*(bitand(d,uint8(4))>0)));
  else
    [height, width] = size(a);
    height_2 = round(height/rendering_scale);
    width_2 = round(width/rendering_scale);
    a = imresize(a, [height_2, width_2]);
    
    image = repmat(a, [1 1 3]);
    d = imresize(diff_stack(:,:,i), [height_2, width_2], 'nearest');
    image(:,:, 1) = ...
      max(image(:,:, 1), uint8(255*(bitand(d,uint8(1))>0)));
    image(:,:, 2) = ...
      max(image(:,:, 2), uint8(255*(bitand(d,uint8(2))>0)));
    image(:,:, 3) = ...
      max(image(:,:, 3), uint8(255*(bitand(d,uint8(4))>0)));
  end
  
  imwrite(image, diff_stack_file_name, 'Compression', 'none', 'WriteMode', 'append');
end
fprintf('done.\n');

return
end
