function [overlap_mask_1, overlap_mask_2] = get_overlap_mask(size_1, size_2, transform)
% [overlap_mask_1, overlap_mask_2] = get_overlap_mask(size_1, size_2, ...
% transform) 
% Compute overlap masks for pair of images given their sizes and the
% transform of the second image with respect to the first one.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  01052009  init. code
%

affine_t = [reshape(transform, [2 3]); 0 0 1];
tform = maketform('affine', (affine_t(1:2, 1:3))');
[image_t, xdata, ydata] = ...
  imtransform(ones(size_2), tform);

minx = min([1, xdata]);
maxx = max([size_1(2), xdata]);
miny = min([1, ydata]);
maxy = max([size_1(1), ydata]);

xdata = round(xdata - minx + 1);
ydata = round(ydata - miny + 1);
maxx = round(maxx - minx + 1);
maxy = round(maxy - miny + 1);

temp_mask = zeros(maxy, maxx);
temp_mask(ydata(1):ydata(2), xdata(1):xdata(2)) = image_t;
xdata_1 = [1, size_1(2)] - minx + 1;
ydata_1 = [1, size_1(1)] - miny + 1;
overlap_mask_1 = temp_mask(ydata_1(1):ydata_1(2), xdata_1(1):xdata_1(2));

affine_t_inv = inv(affine_t);
tform = maketform('affine', (affine_t_inv(1:2, 1:3))');
[image_t, xdata, ydata] = ...
  imtransform(ones(size_1), tform);

minx = min([1, xdata]);
maxx = max([size_2(2), xdata]);
miny = min([1, ydata]);
maxy = max([size_2(1), ydata]);

xdata = round(xdata - minx + 1);
ydata = round(ydata - miny + 1);
maxx = round(maxx - minx + 1);
maxy = round(maxy - miny + 1);

temp_mask = zeros(maxy, maxx);
temp_mask(ydata(1):ydata(2), xdata(1):xdata(2)) = image_t;
xdata_2 = [1, size_2(2)] - minx + 1;
ydata_2 = [1, size_2(1)] - miny + 1;
overlap_mask_2 = temp_mask(ydata_2(1):ydata_2(2), xdata_2(1):xdata_2(2));

return
end
