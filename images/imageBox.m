function [lines mask] = imageBox(im, threshold)
mask = extractImEdge(im, threshold);
[h t r] = hough(mask, 'ThetaResolution', .1);
peaks = houghpeaks(h, 16, 'Threshold', 0);
lines = houghlines(mask, t, r, peaks, 'FillGap', 1000, 'MinLength', 400);
end