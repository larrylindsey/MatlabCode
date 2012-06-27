function [A_NL o_LD] = legendreMat(x, ignore, order)
% [A o] = legendreMat(x, ignore, order)
%  Creates a Vandermonde matrix over the Legendre polynomials.
%
% x - [n k] - input sample locations, for n samples in a k-dimensional space
% order [1] - the maximal order of the power series
%
% A - [n l] - a Vandermonde matrix. For example, given a two dimensional space
%             and order 2, each row in the matrix is
% [P0, P1{x(n, 1)}, P1{x(n, 2)}, P2{x(n,1)}, P1{x(n,1)}*P1{x(n,2)},  P2{x(n,2)}]
%             Where Pk is the normalized kth order Legendre polynomial.
%             Please note the abuse of Matlab notation above.
% o - [l k] - a matrix representing the order of the lth term in the kth
%             dimension. This is the direct output of listOrder(k, order)


if nargin < 3
    order = ignore;
end

[n d] = size(x);

% Lp represents the uni-variate Legendre polynomials for each dimension.
Lp_NOD = zeros(n, order + 1, d);

for ii = 1:d
    Lp_NOD(:,1,ii) = ones(n, 1);
    Lp_NOD(:,2,ii) = x(:,ii);
    for jj = 2:order
        Lp_NOD(:,jj + 1,ii) = nextPolynomial(Lp_NOD(:,jj,ii),...
            Lp_NOD(:,jj - 1, ii), jj, x(:,ii));
    end
end

% o represents the order in each variable for a given monomial. It is indexed
% (l, d), where l indexes the individual monomial, and d is the dimension of the
% given variable. See listOrder() for more details.
o_LD = listOrder(d, order);
l = size(o_LD, 1);

A_NL = zeros(n, l);

for ii = 1:l
	normFactor = prod(sqrt(2 ./ (2 * o_LD(ii,:) + 1)));
    A_NL(:,ii) = Lp_NOD(:, o_LD(ii,1) + 1, 1);
    for jj = 2:d
        A_NL(:,ii) = A_NL(:,ii) .* Lp_NOD(:, o_LD(ii,jj) + 1, jj);
    end
    A_NL(:,ii) = A_NL(:,ii) / normFactor;
end

end

function pn = nextPolynomial(pn_1, pn_2, o, x)

pn = ((2 * o - 1) * x .* pn_1 - (o - 1) * pn_2) / o;

end
