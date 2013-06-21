function eqmap = mergeEqmap(eqmap, bwl3d)

n = size(bwl3d,3);

bwl_min = zeros(n,1);
bwl_max = bwl_min;

bwl_min(end + 1) = inf;

for ii = 1:n
    bwl3d_i = bwl3d(:,:,ii);
    bwl_min(ii) = min(bwl3d_i(bwl3d_i(:) ~= 0));
    bwl_max(ii) = max(bwl3d_i(:));
end

cmap = cell(ceil(n / 2),1);

for ii = 1:floor(numel(bwl_min) / 2)
    sel = and(eqmap(:,1) >= bwl_min((ii - 1) * 2 + 1), ...
        eqmap(:,1) < bwl_min(ii * 2 + 1));
    cmap{ii} = eqmap(sel, 1:2);    
end

if n / 2 ~= floor(n / 2)
    sel = eqmap(:,1) >= bwl_min(end);
    cmap{end} = eqmap(sel, 1:2);
end

fprintf('%d...', numel(cmap));
cmap = cmapConnComp(cmap);

while numel(cmap) > 1    
    cmap = mergeMaps(cmap);
    fprintf('%d...', numel(cmap));
    cmap = cmapConnComp(cmap);    
end
fprintf('\n');

eqmap = cmap{1};

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cmap = cmapConnComp(cmap)
parfor ii = 1:numel(cmap)
    cmap{ii} = graph_cc(cmap{ii});
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mcmap = mergeMaps(cmap)
mcmap = cell(ceil(numel(cmap) / 2), 1);

jj = 1;

for ii = 1:2:(numel(cmap) - 1)
    mcmap{jj} = cat(1, cmap{ii}, cmap{ii + 1});
    jj = jj + 1;
end

if ii + 1 < numel(cmap)
    mcmap{jj} = cmap{end};
end

end
