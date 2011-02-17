function mout = fenceSqDistPost(measure, samples, extra)

L = extra.L;

pkval = zeros(size(samples, 1), 1);

for i_p = 1:numel(pkval)
    pkloc = samples(i_p, 1);
    pkval(i_p) = L(pkloc, 2);
end

pk_extent = max(samples(:,1)) - min(samples(:,1));
fence_select = find(L(:,1) > 0);
fence_extent = max(fence_select) - min(fence_select);

measure = measure + max(L(:,2)) - mean(pkval);

mout = measure * fence_extent * range(pkval) / pk_extent / min(pkval);


end