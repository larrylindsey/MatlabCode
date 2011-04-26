function [ws labelim] = fijiTrainableSegmentationAnnotation(classim)

openstsz = 4;
closestsz = 2;

classim = logical(classim > 0);

classim = imopen(classim, strel('disk', openstsz));
classim = imclose(classim, strel('disk', closestsz));

imd2 = bwdist(classim);

ws = watershed(imd2, 4);

if nargout > 1
    n = max(ws(:));
    cmap = hsv(n);
    cmap = cmap(randperm(n), :);
    cmap = cat(1, zeros(1, 3), cmap);
    
    labelim = applyColorMap(ws, cmap);
end