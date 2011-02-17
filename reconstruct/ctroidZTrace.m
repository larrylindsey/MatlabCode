function ctr = ctroidZTrace(secdoc, cname)

if nargin > 1
    contour = extractContour(secdoc, cname);
else
    contour = secdoc;
end

ctr = zeros(numel(contour), 3);

for ic = 1:numel(contour)
    ctr(ic, :) = contourCentroid(contour(ic));
end


end