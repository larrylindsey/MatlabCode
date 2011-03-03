function [x y] = reconstructDomainBounds(secdoc)

secimtrans = collectReconstructImageTransforms(secdoc);
allContours = [secimtrans.Contour];
allImageDomainTransPts = cat(1, allContours.imageDomainTransPoints);
xymin = min(allImageDomainTransPts, [], 1);
xymax = max(allImageDomainTransPts, [], 1);

x = [xymin(1) xymax(1)];
y = [xymin(2) xymax(2)];

