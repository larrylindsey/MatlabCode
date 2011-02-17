function [c M] = contourCentroid(secdoc, cname)

if nargin > 1
    contour = extractContour(secdoc, cname);
else
    contour = secdoc;
end

section = [contour.section];

z = [contour.z];

m = zeros(size(z));
x = m;
y = m;

for i_sec = 1:numel(section)
    
    currContour = contour(i_sec);
    
    pts = currContour.transPoints;
    n = size(pts, 1);
    pts = pts([1:end 1], :);
    
    m(i_sec) = .5 * sum(pts(1:n,1) .* pts(2:end,2) - ...
        pts(2:end, 1) .* pts(1:n,2));
    
    x(i_sec) = sum( (pts(1:n, 1) + pts(2:end, 1)) .* ...
        (pts(1:n, 1) .* pts(2:end, 2) - pts(2:end, 1) .* pts(1:n, 2))) / ...
        (6 * m(i_sec));

    y(i_sec) = sum( (pts(1:n, 2) + pts(2:end, 2)) .* ...
        (pts(1:n, 1) .* pts(2:end, 2) - pts(2:end, 1) .* pts(1:n, 2))) / ...
        (6 * m(i_sec));
    
end

M = sum(m);

c = [sum(x .* m) sum(y .* m) sum(z .* m)] / M;

end