function cim = applyColorMap(im, cmap, lim)

if size(im, 3) > 1
    im = rgb2gray(im);
end

im = double(im);

if nargin < 3 || isempty(lim)
    im = im - min(im(:));
    im = im / max(im(:));
else
    im = im - lim(1);
    im = im / (lim(2) - lim(1));
    im(im < 0) = 0;
    im(im > 1) = 1;
end

im = floor(im * (size(cmap,1) - 1)) + 1;
im(isnan(im)) = 1;

cimr = zeros(size(im));
cimg = cimr;
cimb = cimr;

cimr(:) = cmap(im,1);
cimg(:) = cmap(im,2);
cimb(:) = cmap(im,3);

cim = cat(3, cimr, cimg, cimb);


end