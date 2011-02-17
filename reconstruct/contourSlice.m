function [slice x y] = contourSlice(secdoc, cname, section, x, y)

if ischar(cname)
    cname = {cname};
end

contour = repmat(struct, 0);

for i_c = 1:numel(cname)
    contour = cat(2, contour, extractContour(secdoc, cname{i_c}));
end

if min(size(x)) == 1
    [x y] = meshgrid(x, y);
end

slice = false(size(x));

if isempty(contour)
    return;
end

csec = [contour.section];

if numel(secdoc) > 1
    csel = csec == section;
else
    csel = 1;
end

contour = contour(csel);

for i_c = 1:numel(contour)
    slice = or(slice, inpolygon(x, y, contour(i_c).transPoints(:,1), ...
        contour(i_c).transPoints(:,2)));
end

end