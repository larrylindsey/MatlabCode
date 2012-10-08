function mask = alphaMask2x2k24(im)

satmask = maskBrightestRegion(im, 20);
im(satmask) = 0;

mask = im == 0;
mask = imclose(mask, strel('disk', 3));

L = bwlabel(mask);
stat = regionprops(L, 'Area');
[~, isel] = max([stat.Area]);

mask = not(L == isel);


mask = imopen(imfill(imclose(mask, strel('disk', 7)), 'holes'), ...
    strel('disk', 5));


mask = not(mask);
end
