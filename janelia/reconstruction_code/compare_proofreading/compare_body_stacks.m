function [diff_stack, diff_stack_blend] = compare_body_stacks(body_stack_a, body_stack_b, image_stack)
% compare_proofread_stacks_from_raveler(...
%   section_image_name_format, planes, proofread_dir_a, proofread_dir_b, ...
%   diff_stack_file_name, rendering_scale, render_original_image)
% Compare two proofreadings of a stack obtained from Raveler.
% The two proofreading sessions must first be exported from Raveler.
%
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  06192009  init. code
% v1  11242009  more general input
%

fprintf('Computing difference between proofread stacks ...\n');
min_overlap_norm_threshold = 0.5;
[~, diff_stack]  = ...
  compare_label_stacks_bp2(body_stack_a, body_stack_b, min_overlap_norm_threshold);
fprintf('done.\n');

if(nargin<=2)
  return;
end

fprintf('Rendering difference stack for viewing ...\n');
diff_stack_blend = cell(size(image_stack,3),1);
for i = 1:size(image_stack,3)
  fprintf('plane: %d\n', i);
  
  image = repmat(image_stack(:,:,i), [1 1 3]);
  image(:,:, 1) = ...
    max(image(:,:, 1), uint8(255*(bitand(diff_stack(:,:,i),uint8(1))>0)));
  image(:,:, 2) = ...
    max(image(:,:, 2), uint8(255*(bitand(diff_stack(:,:,i),uint8(2))>0)));
  image(:,:, 3) = ...
    max(image(:,:, 3), uint8(255*(bitand(diff_stack(:,:,i),uint8(4))>0)));
  diff_stack_blend{i} = image;
end
fprintf('done.\n');

return
end
