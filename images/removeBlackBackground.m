function imout = removeBlackBackground(imin)


if size(imin, 3) > 1
    im = rgb2gray(imin);
else
    im = imin;
end

rsel = max(im, [], 2) > 0;
csel = max(im, [], 1) > 0;

imout = imin(rsel, csel, :);
