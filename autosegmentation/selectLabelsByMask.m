function [maskOut labels weight] = selectLabelsByMask(L, mask, t)

if nargin < 3
    t = .75;
end

maskSize = sum(mask(:));
labels = unique(L(mask));
labels = labels(labels ~= 0);

labelSize = zeros(size(labels));

for il = 1:numel(labels)
    labelSize(il) = sum(L(mask) == labels(il)) / maskSize;
end


[labelSize iOrder] = sort(labelSize, 'descend');
labels = labels(iOrder);

iCutoff = find(cumsum(labelSize) > t, 1, 'first');

labels = labels(1:iCutoff);
weight = labelSize(1:iCutoff);

%labels = cat(1, labels(:), innerLabels(:));
%weight = cat(1, weight(:), innerWeight(:));

maskOut = false(size(L));
for il = 1:numel(labels)
    maskOut = or(maskOut, L == labels(il));
end

