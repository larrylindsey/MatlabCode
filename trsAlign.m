function [opts tvect theta scale] = trsAlign(fpts, tpts)
% function opts = trsAlign(fpts, tpts)
%   Aligns fpts to tpts by translating, rotating, and scaling
%   In other words, after translation, performs a series of linear
%   operations on fpts to minimize the distance to tpts, without affecting
%   the "squareness" of fpts.
%
% fpts [n 2] - the set of points to transform
% tpts [n 2] - the set of points to trasnform onto.
% 
% opts [n 2] - fpts, aligned to tpts.

n = size(fpts, 1);
if size(tpts, 1) ~= n
    error('fpts and tpts must be of the same size');
end

fcom = mean(fpts, 1);
tcom = mean(tpts, 1);
tvect = tcom - fcom;


% translate
opts = fpts - repmat(fcom, [n 1]);
tpts = tpts - repmat(tcom, [n 1]);

% scale
tmag = sqrt(sum(tpts.^2, 2));
fmag = sqrt(sum(opts.^2, 2));
nzsel = tmag > 0 & fmag > 0;
scale = mean(tmag(nzsel) ./ fmag(nzsel));

opts = opts * scale;


% rotate
alpha = sum(sum(opts .* tpts));
beta = sum(opts(:,1) .* tpts(:,2) - opts(:,2) .* tpts(:,1));
theta = asin(sqrt(1 / (alpha * alpha / beta / beta + 1)));
rotmat0 = [cos(theta) sin(theta); -sin(theta) cos(theta)];
rotmat1 = rotmat0';

opts0 = opts * rotmat0;
opts1 = opts * rotmat1;

if mean(sum((opts0 - tpts).^2, 2)) > mean(sum((opts0 + tpts).^2, 2))
    opts0 = -opts0;    
end
if mean(sum((opts1 - tpts).^2, 2)) > mean(sum((opts1 + tpts).^2, 2))
    opts1 = -opts1;
end
if mean(sum((opts0 - tpts).^2, 2)) > mean(sum((opts1 - tpts).^2, 2))
    opts = opts1;
else
    opts = opts0;
end

opts = opts + repmat(tcom, [n 1]);
end