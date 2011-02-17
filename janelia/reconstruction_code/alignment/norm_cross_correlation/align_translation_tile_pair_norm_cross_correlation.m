function transform = align_translation_tile_pair_norm_cross_correlation(...
  image_1, image_2, scale, fold_mask_1, fold_mask_2)
% transform = align_translation_tile_pair_norm_cross_correlation(...
%   image_1, image_2)
% Compute translation based alignment between two images based on
% normalized croos correlation.
% Inputs:
%   image_1, image_2    Images M1xN1, M2xN2 double
%   scale               Scale to be used during alignment for speeding up
% Output:
%   transform .transform   transformation in affine form [1 0 0 1 dx dy]
%             .corr_value  correlation value at maximal peak.
%
% Shiv N. Vitaladevuni
%

if(max(image_1(:))==0 || max(image_2(:))==0)
  transform.transform = [1 0 0 1 0 0];
  transform.corr_value = -1;
  return;
end

image_1_s = imresize(image_1, 1/scale, 'bilinear');
image_2_s = imresize(image_2, 1/scale, 'bilinear');

if(nargin>3)
  fold_mask_1_s = imresize(fold_mask_1, 1/scale, 'nearest');
  fold_mask_2_s = imresize(fold_mask_2, 1/scale, 'nearest');
  image_1_s_n = (fold_mask_1_s>0) .* (image_1_s - mean(mean(image_1_s)));
  image_2_s_n = (fold_mask_2_s>0) .* (image_2_s - mean(mean(image_2_s)));
else
  image_1_s_n = image_1_s - mean(mean(image_1_s));
  image_2_s_n = image_2_s - mean(mean(image_2_s));
end

image_1_s_n_p = zeros(size(image_1_s) + 2*size(image_2_s));
image_1_s_n_p(size(image_2_s_n,1)+1:size(image_2_s_n,1)+size(image_1_s_n,1), ...
  size(image_2_s_n,2)+1:size(image_2_s_n,2)+size(image_1_s_n,2)) = image_1_s_n;
  
try
	cc = normxcorr2(image_2_s_n, image_1_s_n_p);
catch
  transform.transform = [1 0 0 1 0 0];
  transform.corr_value = -1;
  return;
end	
cc = cc(size(image_2_s_n,1)+1:end, size(image_2_s_n,2)+1:end);

[max_val, max_id] = max(cc(:));
[peak_y, peak_x] = ind2sub(size(cc), max_id);
c_x_s = peak_x - size(image_2_s_n, 2);
c_y_s = peak_y - size(image_2_s_n, 1);
c_x = c_x_s * scale;
c_y = c_y_s * scale;

if(c_x>=size(image_1,2) || c_y>=size(image_1,1) || ...
    c_x<=-size(image_2,2) || c_y<=-size(image_2,1))
  max_val = -1;
end

transform.transform = [1 0 0 1 c_x, c_y];
transform.corr_value = max_val;

return
end
