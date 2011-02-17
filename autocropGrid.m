function cropsel = autocropGrid(edgeMap)

maxDim = 4096;

edgePeaks = findPeaks2(edgeMap);
edgePeakValSort = sort(edgeMap(edgePeaks));

maxVal = edgePeakValSort(end);

valSel = edgePeakValSort >= .01 * maxVal & edgePeakValSort <= .99 * maxVal;
edgePeakValSort = edgePeakValSort(valSel);

minVal = edgePeakValSort(round(end * .25));
maxVal = edgePeakValSort(round(end * .75));

edgePeaks = edgePeaks & edgeMap >= minVal & edgeMap <= maxVal;
[r c] = find(edgePeaks);

resize_fact = max(1, max(size(edgeMap) / maxDim));
r = r / resize_fact;
c = c / resize_fact;
sz = round(size(edgeMap) / resize_fact);

k = convhull(r, c);

[cc rr] = meshgrid(1:sz(1), 1:sz(2));

cropsel = inpolygon(rr, cc, r(k), c(k));

cropsel = imresize(cropsel, mean(size(edgeMap) ./ sz), 'nearest');



% d = bwdist(edgePeaks);
% 
% cropsel = d < median(d(:));
