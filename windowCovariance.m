function imcov = windowCovariance(im1, im2, n)

if numel(n) == 1
    k = ones(n) / n / n;
else
    k = n;
end

E1 = imfilter(im1, k);
E2 = imfilter(im2, k);
E12 = imfilter(im1.*im2, k);

imcov = E12 - E1.*E2;
