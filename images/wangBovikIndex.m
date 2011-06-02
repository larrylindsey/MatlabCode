function [Q Qcat] = wangBovikIndex(im1, im2, n, p)

if nargin < 4
    p = false;
end

im1 = im2single(im1);
im2 = im2single(im2);

if size(im1, 3) > 1
    im1 = rgb2gray(im1);
end

if size(im2, 3) > 1
    im2 = rgb2gray(im2);
end

if p
    parIm = cat(3, im1, im2, im1.*im2, im1, im2);
    parFun = {@meanFun @meanFun @meanFun @varFun @varFun};
    parImout = zeros(size(parIm));
    
    parfor ii = 1:5
        parImout(:,:,ii) = slidingImageWindow(parIm(:,:,ii), ...
            [n n], parFun{ii});
    end
    
    clear parIm;
    mu1 = parImout(:,:,1);
    mu2 = parImout(:,:,2);
    mu12 = parImout(:,:,3);
    dev1 = parImout(:,:,4);
    dev2 = parImout(:,:,5);
    clear parImout;
else
    mu1 = slidingImageWindow(im1, [n n], @meanFun);
    mu2 = slidingImageWindow(im2, [n n], @meanFun);
    mu12 = slidingImageWindow(im1.*im2, [n n], @meanFun);
    
    dev1 = slidingImageWindow(im1, [n n], @varFun);
    dev2 = slidingImageWindow(im2, [n n], @varFun);
end

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

function m = meanFun(im)
m = mean(im(:));
end

function v = varFun(im)
v = mean(im(:).^2) - mean(im(:)).^2;
end