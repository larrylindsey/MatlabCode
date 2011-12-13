function gridDistortionAffineFigs(prefix, grid_model, rc_found, rc_grid,...
    norm)

fmt = {'png', 'fig'};
order = [1 2 3 4 5 6 7];
trstr = extractedTransform(rc_found, rc_grid, grid_model, order);


clf;
imageGridDifference(rc_grid, rc_found, trstr.affine(1).rc_grid, norm / 2);
savefig(prefix, 'found', fmt);

for ii = 1:numel(order)
    clf;
    imageDistortionMap(trstr.affine(ii).tr, [100 50], [-1 1], norm / 2);
    savefig(prefix, sprintf('dmap_o%d', order(ii)), fmt);
end

end

function savefig(prefix, f, fmt)

for ii = 1:numel(fmt)
    saveas(1, sprintf('%s_%s.%s', prefix, f, fmt{ii}));
end

end