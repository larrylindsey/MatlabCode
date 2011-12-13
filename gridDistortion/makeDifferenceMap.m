function [dfmap_mag dfmap_color rc_grid] =...
    makeDifferenceMap(rc_grid, rc1, rc2, n)

if nargin < 4
    n = [];
end

if size(rc_grid, 2) > size(rc_grid, 1)
    rc_grid = rc_grid.';
end

rc_grid(:,1) = rc_grid(:,1) + 1 - min(rc_grid(:,1));
rc_grid(:,2) = rc_grid(:,2) + 1 - min(rc_grid(:,2));

sz = max(rc_grid, [], 1);
dfmap_mag = nan(sz);
dfmap_color = repmat(dfmap_mag, [1 1 3]);

rcDifference = rc1 - rc2;

dfmap_mag(sub2ind(sz, rc_grid(:,1), rc_grid(:,2))) = ...
    sqrt(sum(rcDifference.^2, 2));

imv = vectorMapToImage(rcDifference(:,2), rcDifference(:,1), n);
dfmap_color(sub2ind(size(dfmap_color), rc_grid(:,1), rc_grid(:,2), ...
    ones(size(rc_grid(:,1))))) = imv(:,1);
dfmap_color(sub2ind(size(dfmap_color), rc_grid(:,1), rc_grid(:,2), ...
    2*ones(size(rc_grid(:,1))))) = imv(:,2);
dfmap_color(sub2ind(size(dfmap_color), rc_grid(:,1), rc_grid(:,2), ...
    3*ones(size(rc_grid(:,1))))) = imv(:,3);