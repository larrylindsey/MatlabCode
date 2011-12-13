function imageGridDifference(rc_grid, rc1, rc2, n, varargin)
tn = [];
[t varargin] = getargs({'threshold'}, varargin);
t = t.threshold;
if ~isempty(t)
    tn = t / n;
end


[dfmag dfcolor rc_grid] = makeDifferenceMap(rc_grid, rc1, rc2, tn);

figure;
imagesc(dfmag * n); axis image; axis xy;
if ~isempty(t)
    caxis([0 t]);
end
colormap(jet(256));
colorbar;

figure;
imagesc(dfcolor); axis image; axis xy;
hold on;
quiver(rc_grid(:,2), rc_grid(:,1), ...
    rc1(:,2) - rc2(:,2), rc1(:,1) - rc2(:,1), varargin{:})