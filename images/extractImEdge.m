function edgeIm = extractImEdge(im, threshold)

if nargin < 2
    threshold = 0;
end

mask = imclose(im, strel('disk', 32)) <= threshold;

maskOutter = imerode(mask, strel('disk', 2));

edgeIm = mask - maskOutter;