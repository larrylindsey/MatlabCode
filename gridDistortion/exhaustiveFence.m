function [sel bestMod measureOut vout] = exhaustiveFence(uni, ex)
n = size(uni, 1);
measureOut = nan(n-1, n);
vout = measureOut;

for ii=1:(n-1)
for jj = (ii+1):n
    model = fenceModel(uni([ii jj], :), []);
    [measureOut(ii,jj) vout(ii, jj)] = exhaustiveFenceMetric(model, ex);    
end
end

mbest = min(measureOut(:));
[pp qq] = find(measureOut == mbest);
bestMod = zeros(numel(pp), 2);

for ii = 1:numel(pp)
    bestMod(ii, :) = fenceModel(uni([pp(ii) qq(ii)], :));
end

bestMod = unique(bestMod, 'rows');

labelMaster = 1:max(ex(:,2));
sel = zeros(size(bestMod, 1), numel(labelMaster));
for i_m = 1:size(bestMod, 1)
    v = bestMod(i_m, 1):bestMod(i_m, 2):size(ex, 1);
    L = ex(v,2);
    for i_l = 1:numel(labelMaster)
        sel(i_m, i_l) = any(L == labelMaster(i_l));
    end
end

sel = unique(sel, 'rows');

if size(sel, 1) > 1
    warning('Fence:discriminate',['Could not discriminate between best'...
        ' models based on selection.  Attempting to do it by quantity']);
    sval = sum(sel, 2);
    smin = min(sval);
    i_smin = logical(sval == smin);
    sel = sel(i_smin, :);
    if size(sel, 1) > 1
        warning('Fence:discrimitate', ['Could not discriminate by ' ...
            'quantity.  Attempting to do it by spacing']);
        mmax = max(bestMod(:, 2));
        i_mmax = logical(bestMod(:, 2) == mmax);
        sel = sel(i_mmax, :);
        if size(sel, 1) > 1
            warning('Fence:discriminate', ['Could not discriminate by '...
                'spacing.  Just picking the first one.']);
            sel = sel(1,:);
        end
    end
end

