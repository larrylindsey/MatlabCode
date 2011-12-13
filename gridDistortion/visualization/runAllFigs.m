function runAllFigs(matfiles, order, figroot)

if nargin < 2
    order = 5;
end

for ii = 1:numel(matfiles)
    dstr(ii) = load(matfiles{ii}, 'rc_found', 'rc_grid', 'scale',...
        'grid_model', 'fname');
    trstr(ii) = extractedTransform(dstr(ii).rc_found, dstr(ii).rc_grid, ...
        dstr(ii).grid_model, order);
    magmap = makeDifferenceMap(dstr(ii).rc_grid, dstr(ii).rc_found, ...
        trstr(ii).affine(1).rc_grid); % mean(dstr(ii).scale / 2)
    
    affMag(ii) = max(magmap(:) * mean(dstr(ii).scale / 2));
    
    magmap = makeDifferenceMap(dstr(ii).rc_grid, dstr(ii).rc_found, ...
        trstr(ii).square(1).rc_grid); % mean(dstr(ii).scale / 2)
    
    sqMag(ii) = max(magmap(:) * mean(dstr(ii).scale / 2));
end

unix(['mkdir -p ' figroot]);

affMag = ceil(max(affMag));
sqMag = ceil(max(sqMag));

close all;

for ii = 1:numel(matfiles)
    
    idot = find(dstr(ii).fname == '.', 1, 'first');
    name = upper(dstr(ii).fname(1:(idot - 1)));
    
    imageGridDifference(dstr(ii).rc_grid, dstr(ii).rc_found, ...
        trstr(ii).square(1).rc_grid, mean(dstr(ii).scale / 2), ...
        'threshold', sqMag);
    imageGridDifference(dstr(ii).rc_grid, dstr(ii).rc_found, ...
        trstr(ii).affine(1).rc_grid, mean(dstr(ii).scale / 2), ...
        'threshold', affMag);
    
    writeGridFig(name, figroot, 'Orthogonal', 1, 2);
    writeGridFig(name, figroot, 'Affine', 3, 4);
    close all;
    
    for jj = 1:numel(order)
        imageDistortionMap(trstr(ii).affine(jj).tr, ...
            [512 ceil(sqrt(size(dstr(ii).rc_found,1)))], ...
            [-1 1], mean(dstr(ii).scale) / 2);
        imageDistortionMap(trstr(ii).affine(jj).tr, ...
            [512 ceil(sqrt(size(dstr(ii).rc_found,1)))], ...
            [-1 1], mean(dstr(ii).scale) / 2, 'rc', dstr(ii).rc_found, ...
            'threshold', affMag);
        writeTransformFig(name, figroot, order(jj), true, 1, 2);
        writeTransformFig(name, figroot, order(jj), false, 3, 4);
        close all;
    end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function formatTitle(fig, tcap)
figure(fig);
title(tcap, 'FontSize', 16, 'Interpreter', 'none');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dosave(fig, prefix)
fmt = {'png', 'fig'};

for ff = 1:numel(fmt)
    saveas(fig, [prefix '.' fmt{ff}], fmt{ff});
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeGridFig(name, figroot, adj, f1, f2)
tcap = sprintf([adj ' Grid Crossing\nDisplacement Magnitude for ' name]);
formatTitle(f1, tcap);
tcap = sprintf([adj ' Grid Crossing\nDisplacement Vector Map for ' name]);
formatTitle(f2, tcap);

dosave(f1, sprintf('%s/%s_%s_grid_mag', figroot, name, adj));
dosave(f2, sprintf('%s/%s_%s_grid_vmap', figroot, name, adj));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeTransformFig(name, figroot, o, iscropped, f1, f2)
tcap = sprintf('Transform Approximation Magnitude\nAt Order %d for %s',...
    o, name);
formatTitle(f1, tcap);
tcap = sprintf('Transform Approximation Vector Map\nAt Order %d for %s',...
    o, name);
formatTitle(f2, tcap);

if iscropped
    cstr = '_cropped';
else
    cstr = '';
end

dosave(f1, sprintf('%s/%s_trans_mag_o%d%s', figroot, name, o, cstr));
dosave(f2, sprintf('%s/%s_trans_vmap_o%d%s', figroot, name, o, cstr));
end