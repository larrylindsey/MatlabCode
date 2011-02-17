function fold_mask = compute_fold_mask(image, config)
% fold_mask = compute_fold_mask(image)
% Compute fold mask to assist segmentation and alignment
%
% Inputs:
%   image       uint8 matrix of the image
% Output
%   fold_mask   uint8 mask with fold pixels marked 0, and connected
%                 components labelled with integers.
%
% Shiv N. Vitaladevuni
% 

if(nargin>1)
  min_fold_area = config.min_fold_area;
else
  min_fold_area = 1000;
end

fold_mask_init = get_fold_mask(image);

fold_mask_a = bwareaopen(fold_mask_init, min_fold_area);

fold_mask_a_f = imfill(fold_mask_a, 'holes');

mf_element = strel('disk', 10);
fold_mask_a_f_de = imerode(imdilate(fold_mask_a_f, mf_element), mf_element);

fold_mask = bwlabel(fold_mask_a_f_de==0);

return
end
