function y2 = doRadialSpiralTransform(T, x2, o, ~, pfun)

fun = pfun.fun;

if size(x2, 2) ~= 2
    error('radialAndSpiralMat supports only 2 dimensions, rather than %d', ...
        size(x2, 2));
end

A2 = radialAndSpiralMat(x2, [], o, fun);

u = A2(:,:,1) * T(:,1);
v = A2(:,:,2) * T(:,2);

y2 = cat(2, u, v);

end
