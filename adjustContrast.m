function imOut = adjustContrast(imIn, ranges)
% function imOut = adjustContrast(imIn, [ranges])
% Scales the values of a given image so that the output uses the full
% dynamic range.
%
% imIn - the image whose contrast is to be scaled.
% ranges - controls how minimum and maximum values are sampled.  This is
%          a simple operation in practice, but difficult to explain
%          briefly - I suggest you read the code.
%          
% imOut - the contrast-scaled output image.

if nargin < 2
    ranges = [0 .95];
else
    ranges(1) = max(ranges(1), 0);
    ranges(2) = min(ranges(2), 1);
end

imOut = double(imIn);
imGray = mean(imOut, 3);
n = numel(imGray);
sortedIm = sort(imGray(:));

minIndex = round(n * ranges(1));
minIndex = max(1, minIndex);
minVal = sortedIm(minIndex);
minVal = max(0, minVal);

imOut = imOut - minVal;
imOut(logical(imOut < 0)) = 0;

imGray = mean(imOut, 3);
sortedIm = sort(imGray(:));

maxIndex = round(n * ranges(2));
maxIndex = min(maxIndex, n);
maxVal = sortedIm(maxIndex);
maxVal = min(maxVal / ranges(2), sortedIm(end));

imOut = imOut / maxVal;
imOut(logical(imOut > 1)) = 1;
