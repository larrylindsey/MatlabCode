function [A o] = rationalMat(x, y, order, vm)

if nargin < 1
    A = @doRationalTransform;
    o = @dimWiseLeastSquaresFit;
    return;
end

if nargin < 4
    vm = @taylorMat;
elseif isstruct(vm)
    vm = vm.fun;
end

if numel(order) == 1
    sprintf(['Singleton order.', ...
        ' Assuming equal nominator/denominator order = %d\n'], order);
    order = [order order];
end

yd = size(y, 2);

[P o_n] = vm(x, [], order(1));
[Q o_d] = vm(x, [], order(2));
Q = -Q;

Q(:,1) = [];

Qd = Q .* repmat(y(:,1), [1 size(Q, 2)]);
A = cat(2, P, Qd);
A = repmat(A, [1 1 yd]);

for ii = 2:yd
    Qd = Q .* repmat(y(:,ii), [1 size(Q, 2)]);
    A(:,:,ii) = cat(2, P, Qd);
end

o = cat(1, o_n, -o_d);
