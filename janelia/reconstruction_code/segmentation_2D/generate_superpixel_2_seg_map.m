function generate_superpixel_2_seg_map(config)
% generate_superpixel_2_seg_map(config)
% generate superpixel_2_seg_map from superpixel_to_seg_label obtained from
% segmentation
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

if(isfield(config.stack, 'image_structure') && ~isempty(config.stack.image_structure))
  generate_superpixel_2_seg_map_multi_tile(config);
else
  generate_superpixel_2_seg_map_uni_tile(config);
end

return
end

