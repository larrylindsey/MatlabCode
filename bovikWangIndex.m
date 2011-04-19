function [Q Qcat] = bovikWangIndex(im1, im2, n)

im1 = im2single(im1);
im2 = im2single(im2);

if size(im1, 3) > 1
    im1 = rgb2gray(im1);
end

if size(im2, 3) > 1
    im2 = rgb2gray(im2);
end

k = ones(n);
k = k / sum(k(:));

mu1 = imfilter(im1, k);
mu2 = imfilter(im2, k);
im1sq = imfilter(im1.*im1, k);
im2sq = imfilter(im2.*im1, k);
mu12 = imfilter(im1.*im2, k);

dev1 = im1sq - mu1.^2;
dev2 = im2sq - mu2.^2;


sigma1 = sqrt(dev1);
sigma2 = sqrt(dev2);
sigma12 = mu12 - mu1.*mu2;

Q1 = sigma12 ./ (sigma1 .* sigma2);
Q2 = 2 * mu1 .* mu2 ./ (mu1.^2 + mu2.^2);
Q3 = 2 * sigma1 .* sigma2 ./ (dev1 + dev2);

Q = Q1 .* Q2 .* Q3;

if nargout > 1
    Qcat = cat(3, Q1, Q2, Q3);
end

end