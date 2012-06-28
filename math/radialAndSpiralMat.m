function [A o] = radialAndSpiralMat(x2, ~, order, fun)

if nargin < 1
    A = @doRadialSpiralTransform;
    o = @dimWiseLeastSquaresFit;
    return;
end

if nargin < 4
    fun = @taylorMat;
elseif isstruct(fun)
    fun = fun.fun;
end

if size(x2, 2) ~= 2
    error('radialAndSpiralMat supports only 2 dimensions, rather than %d', ...
        size(x2, 2));
end

x = x2(:,1);
y = x2(:,2);

r = sqrt(sum(x2.^2, 2));

[P op] = fun(r, order);
[S os] = fun(r, 2);

l = size(P,2);

P1 = repmat(x, [1 l]) .* P;
P2 = repmat(y, [1 l]) .* P;
S1 = repmat(-y, [1 3]) .* S;
S2 = repmat(x, [1 3]) .* S;


A1 = cat(2, P1, S1);
A2 = cat(2, P2, S2);

A = cat(3, A1, A2);
o = cat(1, op, os) + 1;

end
