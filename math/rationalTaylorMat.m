function [A o] = rationalTaylorMat(x, y, order, vm)

if nargin < 4
    vm = @taylorMat;
end

if numel(order) == 1
    sprintf(['Singleton order.', ...
        ' Assuming equal nominator/denominator order = %d\n'], order);
    order = [order order];
end

[P o_n] = vm(x, [], order(1));
[Q o_d] = vm(x, [], order(2));
Q = -Q;

Q(:,1) = [];

Q = Q .* repmat(y(:), [1 size(Q, 2)]);

A = cat(2, P, Q);
o = cat(1, o_n, -o_d);
