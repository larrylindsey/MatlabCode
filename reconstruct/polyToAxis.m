function [a pt] = polyToAxis(c, pts, doconvhull)

if nargin < 3
    doconvhull = false;
end

pts = pts - repmat(c, [size(pts, 1), 1]);

if doconvhull
    k = convhull(pts(:,1), pts(:,2));
else
    k = 1:size(pts, 1);
end

%kpts = pts(k,:);
n = numel(k);
maxMag = max(sqrt(sum(pts.^2, 2))) * 10;

edgeIndices = zeros(2, n);
edgeIndices(1,:) = 1:n;
edgeIndices(2,:) = [2:n 1];
ka = zeros(1, n);
outpt = zeros(n, 2);
kpts = pts(k,:);

parfor ip = 1:n

    currPt = kpts(ip,:);
    rayEnd = -currPt * maxMag;
    axisRay = cat(1, currPt, rayEnd);    
    width = 0;
    for il = 1:n
        polyseg = cat(1, kpts(edgeIndices(1, il), :), ...
            kpts(edgeIndices(2,il), :));
        [intpt doCross] = intersection(axisRay', polyseg');
        intpt = intpt';
        if doCross
            testWidth = sqrt(sum((currPt - intpt).^2, 2));
            if testWidth > width
                width = testWidth;
                outpt(ip, :) = intpt + c;
            end
        end
    end
    ka(ip) = width;
end

a = nan(1, size(pts, 1));
a(k) = ka;

pt = nan(size(pts, 1), 2);
pt(k,:) = outpt;

end

