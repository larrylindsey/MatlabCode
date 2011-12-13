function imout = mergeAlignedImages(im1, im2, x1, y1, x2, y2)

x1 = round(x1);
y1 = round(y1);

xmin = min(cat(1, x1(:), x2(:)));
xmax = max(cat(1, x1(:), x2(:)));
ymin = min(cat(1, y1(:), y2(:)));
ymax = max(cat(1, y1(:), y2(:)));

xmax = xmax + 1 - xmin;
ymax = ymax + 1 - ymin;

imout = zeros(ymax, xmax);

imout((x1(1):x1(2)) + 1 - xmin, (y1(1):y1(2)) + 1 - ymin) = im1(:,:,1);
imout((x2(1):x2(2)) + 1 - xmin, (y2(1):y2(2)) + 1 - ymin) = im2(:,:,1);


end