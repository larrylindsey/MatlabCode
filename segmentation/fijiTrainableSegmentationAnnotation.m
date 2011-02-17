function [ws labelim] = fijiTrainableSegmentationAnnotation(classim, w)

openstsz = 4;
closestsz = 2;

classim = logical(classim > 0);
classim2 = classim;

classim = imopen(classim, strel('disk', openstsz));
classim = imclose(classim, strel('disk', closestsz));

%keyboard

imd2 = bwdist(classim);

ws = watershed(imd2, 4);

if nargin > 1
    ws_sel = ws == 0;
    ws_sel = imdilate(ws_sel, strel('disk', w));
    ws(ws_sel) = 0;
end

if nargout > 1
    
    labelim = zeros([size(classim) 3]);
    L = max(ws(:));
    
    cmap = hsv(L);
    
    cmap = cmap(randperm(L), :);
    cmap = cat(1, zeros(1, 3), cmap);
    
    for il = 1:L
        [r c] = find(ws == il);
        for ic = 1:3
            labelim(sub2ind(size(labelim), r, c, ic * ones(size(r)))) = ...
                cmap(il, ic);
        end
    end
end