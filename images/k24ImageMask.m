function [cnt m b] = k24ImageMask(file, th)


im = imread(file);

if nargin < 2
    im0sort = sort(im(:));
    th = im0sort(round(end / 4));
    clear im0sort;
end

nzmask = im > th;
nzsel = find(nzmask);

if isempty(nzsel)
    fprintf('Found empty image\n');
    cnt = 0;
    m = 0;
    b = 0;
    return;
end

[imsort isort] = sort(im(nzmask));

n = numel(imsort);
x1 = round(n * .1);
x2 = round(n * .9);
fitval1 = imsort(x1);
fitval2 = imsort(x2);

m = double(fitval2 - fitval1) / double(x2 - x1);
b = -m*double(x1) + double(fitval1);

fitline = m * (1:numel(imsort)) + b;
vv = var(double(imsort) - fitline');

sel = and(imsort' < fitline + sqrt(vv), imsort' > fitline - sqrt(vv));

thmin = min(im(nzsel(isort(sel))));
thmax = max(im(nzsel(isort(sel))));

%mask = false(size(im));
%mask(nzsel(isort(sel))) = true;
mask = and(im >= thmin, im <= thmax);

mask = imopen(imclose(mask, strel('disk', 8)), strel('disk', 8));

mask = not(bwareaopen(not(mask), round(numel(mask) / 16)));
mask = bwareaopen(mask, round(numel(mask) / 16));

cnt = sum(mask(:));

imwrite(mask, [file '.mask.png']);
imwrite(compositeImage(im2double(im), mask), [file '.compositeMask.png']);

end