function [r rxy] = calculateAlignmentR(n, pts1, pts2, im1, im2, s, order)

if nargin < 7
    order = 1;
    if nargin < 6
        s = 2;
    end
end

if size(im1, 3) > 1
    im1 = rgb2gray(im1);
end

if size(im2, 3) > 1
    im2 = rgb2gray(im2);
end

im1 = im2single(im1);
im2 = im2single(im2);

tr = regressionTransform(pts1, pts2, order);
tr.data.u = [1 size(im1,2)];
tr.data.v = [1 size(im1,1)];

disp('Applying transform to image data');

im1tr = applyTransformImage(im1, tr, [1 size(im2, 2)], [1 size(im2, 1)]);
im1trmask = applyTransformImage(true(size(im1)), ...
    tr, [1 size(im2, 2)], [1 size(im2, 1)], 'nearest');

f = fspecial('gaussian', ceil(s * 5), s);

disp('Blurring images');
im1tr = imfilter(im1tr, f);
im2 = imfilter(im2, f);

c = clock;
str = sprintf('Rinput_%d_%d_%d_%d_%g', c(2), c(3), c(4), c(5), c(6));

imwrite(im1tr, [str '_01tr.jpg']);
imwrite(im2, [str '_02.jpg']);

im1tr(not(im1trmask)) = nan;
r = corrcoef(im1tr(im1trmask), im2(im1trmask));

disp('Calculating R map');
tic;
if ~isempty(n)
    rxy = rmap(im1tr, im2, n);
end
toc;

end