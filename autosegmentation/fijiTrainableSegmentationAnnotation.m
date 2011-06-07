function [ws labelim] = fijiTrainableSegmentationAnnotation(classim,...
    openstsz, closestsz)

if nargin < 3 || isempty(closestsz)
    closestsz = 2;
end
if nargin < 2 || isempty(openstsz)
    openstsz = 4;
end


classim = logical(classim > 0);

if openstsz > 0
    classim = imopen(classim, strel('disk', openstsz));
end

if closestsz > 0
    classim = imclose(classim, strel('disk', closestsz));
end

imd2 = bwdist(classim);

ws = watershed(imd2, 4);

if nargout > 1
    n = max(ws(:));
    cmap = hsv(n);
    cmap = cmap(randperm(n), :);
    cmap = cat(1, zeros(1, 3), cmap);
    
    labelim = applyColorMap(ws, cmap);
end