function [a idx] = contourAreaBySection(contour, secdoc)

if ischar(contour)
    contour = extractContour(secdoc, contour);
end

idx = unique([contour.section]);
a = zeros(size(idx));

for ii = 1:numel(contour)
    area = polyarea(contour(ii).transPoints(:,1), contour(ii).transPoints(:,2));
    sel = find(contour(ii).section == idx);
    a(sel) = a(sel) + area;
end
