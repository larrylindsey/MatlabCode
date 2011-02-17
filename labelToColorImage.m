function labelim = labelToColorImage(l)

labelim = zeros([size(l) 3]);
L = max(l(:));

cmap = hsv(L);

cmap = cmap(randperm(L), :);
cmap = cat(1, zeros(1, 3), cmap);

for il = 1:L
    [r c] = find(l == il);
    for ic = 1:3
        labelim(sub2ind(size(labelim), r, c, ic * ones(size(r)))) = ...
            cmap(il, ic);
    end
end
