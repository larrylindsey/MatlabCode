function [closest] = expectedGrid(r, c, rTest, cTest)

r = r(:);
c = c(:);
rTest = rTest(:);
cTest = cTest(:);

vect1 = cat(2, r, c);
vect2 = cat(2, rTest, cTest);

dd = dist2(vect1, vect2);

[junk iClosest] = min(dd, [], 2);

closest = vect2(iClosest, :);
end