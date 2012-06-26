function xys = doRadialAndSpiralTransform(xy, tr)

if isstruct(tr)
    coef = tr.T;
else
    coef = tr;
end

order = size(coef, 1) - 4;

[A1 A2] = radialAndSpiralMat(xy(:,1), xy(:,2), order);

xs = A1 * coef(:,1);
ys = A2 * coef(:,2);

xys = cat(2, xs, ys);
