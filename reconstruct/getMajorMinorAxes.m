function [pt d apt] = getMajorMinorAxes(c, pts, convex)
% function [pt d apt] = getMajorMinorAxes(c, pts, [convex])
% c - centroid
% pts - contour points.
% convex - true for convex axes, false otherwise.
%
% pt - points corresponding to major/minor axes.
%      pt(1,:) - major axis
%      pt(2,:) - minor axis
% d - diameters corresponding to major/minor axes.
% apt - points corresponding to the opposite ends of the major and minor
%       axes.
%       

if nargin < 3
    convex = false;
end

[a apts] = polyToAxis(c, pts, convex);
[d imax] = max(a);
[d(2) imin] = min(a);

if min(a) == 0
    keyboard;
end

pt = pts(imax,:);
pt(2,:) = pts(imin,:);

apt = apts([imax imin], :);

