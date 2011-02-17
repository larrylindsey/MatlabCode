function ptsOut = mergeContour(pts1, pts2)

dmap = dist2(pts1, pts2);

dmin = min(dmap(:));

[i1 i2] = find(dmap == dmin);

i1 = i1(1);
i2 = i2(1);

n1 = size(pts1, 1);
n2 = size(pts2, 1);

sel1 = mod((1:n1) + i1 - 1 , n1);
sel1(sel1 == 0) = n1;

sel2 = mod((1:n2) + i2, n2);
sel2(sel2 == 0) = n2;

ptsOut = cat(1, pts1(sel1, :), pts2(sel2, :));

end