function partLabel = fillNN(partLabel, segMask, cutMask)
idL = unique(partLabel(cutMask));

idL(idL == 0) = [];

[rr cc] = find(segMask);
minr = min(rr); maxr = max(rr); minc = min(cc); maxc = max(cc);
segMaskCrop = segMask(minr:maxr, minc:maxc);
%cutMaskCrop = cutMask(min(rr):max(rr), min(cc):max(c));
partLabelCrop = partLabel(minr:maxr, minc:maxc);

bwDistStack = zeros([size(segMaskCrop) numel(idL)]);

for ii = 1:numel(idL)
    bwDistStack(:,:,ii) = bwdist(partLabelCrop == idL(ii));
end

[junk iNN] = min(bwDistStack, [], 3);

nnLabel = idL(iNN);

partLabel(segMask) = nnLabel(segMaskCrop);

end