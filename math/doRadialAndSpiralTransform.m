function xys = doRadialAndSpiralTransform(xy, coef)

order = size(coef, 1) - 4;

[A1 A2] = radialAndSpiralMat(xy, order);

xs = A1 * coef(:,1);
ys = A2 * coef(:,2);

xys = cat(2, xs, ys);
