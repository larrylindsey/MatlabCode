function [A o] = taylorMatND(x, ignore, order)
% [A o] = taylorMatND(x, [~,] order)
%
% Creates a Vandermonde matrix over the power series.
%
% x - [n k] - input sample locations, for n samples in a k-dimensional space
% order [1] - the maximal order of the power series
%
% A - [n l] - a Vandermonde matrix. For example, given a two dimensional space
%             and order 2, each row in the matrix is
%             [1 x(n, 1) x(n, 2) x(n,1)^2 x(n,1)*x(n,2) x(n,2)^2]
% o - [l k] - a matrix representing the order of the lth term in the kth
%             dimension. This is the direct output of listOrder(k, order)

if nargin < 3
    order = ignore;
end

[n k] = size(x);

%o = (order + 1) * (order + 2) / 2;
o = listOrder(k, order);
l = size(o, 1);

A = zeros(n, l);

for ii = 1:l
    A(:,ii) = x(:,1).^o(ii,1);
    for kk = 2:k
        A(:,ii) = A(:,ii).*(x(:,kk).^o(ii,kk));
    end
end
