function plotDMAP(dmap, clim, cmap)

figure;

imagesc(dmap); axis image;
caxis(clim);
colormap(cmap);
cbh = cblabel('Mismatch (pixel)', 'FontSize', 14, 'FontName', 'Arial', ...
    'FontWeight', 'Bold');

set(cbh, 'FontSize', 12);
set(cbh, 'FontName', 'Arial');
set(gca, 'XTick', []);
set(gca, 'YTick', []);


