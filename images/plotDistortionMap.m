function fout = plotDistortionMap(rc_found, rc_grid, sc, order)

% Right now, only handle square images.

sc = sc(1);
rc_fit_sq = trsAlign(rc_grid, rc_found);
rc_fit_aff = affineAlign(rc_grid, rc_found);

data.u = [-1 1];
data.v = [-1 1];
data.n = 100;

ph(1) = doPlot(rc_fit_sq, rc_found, sc, order, data, 'Square');
ph(2) = doPlot(rc_fit_aff, rc_found, sc, order, data, 'Affine');



if nargout > 0
    fout = ph;
end

end

function ph = doPlot(rcfit, rc, sc, order, data, titlestr)
ph.fh = figure;
tr = regressionTransform(rc, rcfit, order, @legendreMat, data);
dmap = makeDistortionMap(tr, 1024, [-1 1]);
ph.imh = imagesc(dmap * sc, 'XData', [1 sc], 'YData', [1 sc]);
axis image;
hold on;
colormap(jet(512));

[r c] = meshgrid(linspace(-1, 1, 15) * .95, linspace(-1, 1, 15) * .95);

rcmap = cat(2, r(:), c(:));
rctrmap = doTransform(rcmap, tr);
rctrd = rctrmap - rcmap;

disp_rcmap = sc * (rcmap + 1) / 2;

ph.qh = quiver(disp_rcmap(:,1), disp_rcmap(:,2), rctrd(:,1), rctrd(:,2),...
    'Color', [1 1 1]);

ph.th = title(titlestr);

ph.cx = caxis;

set(ph.fh, 'UserData', tr);

colorbar;
end