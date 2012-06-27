function [A_NL o_LD] = legendreMat(x, order)

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


% function [Aout x_order y_order] = legendreMat(x, y, order)
% 
% if nargin < 3
%     order = 3;
% end
% 
% x = x(:);
% y = y(:);
% 
% t = numel(x);
% o = (order + 1) * (order + 2) / 2;
% s = 0:order;
% 
% if order < 0
%     Aout = zeros(t, 0);
%     return;
% elseif order == 0
%     Aout = ones(t, 1) * .25;
%     return;
% end
% 
% Aout = zeros(t, o);
% 
% y_order = zeros(1, o);
% x_order = zeros(1, o);
% x_order(1:3) = [0 1 0];
% y_order(1:3) = [0 0 1];
% 
% Aout(:,1) = ones(t, 1);
% Aout(:,2) = x;
% Aout(:,3) = y;
% 
% i_xbasis = s .* (s + 1) / 2 + 1;
% i_ybasis = (s + 2) .* (s + 1) / 2;
% 
% 
% for i_s = 3:numel(s)
%     n = s(i_s);
%     
%     %x functions
%     ix_pn_1 = i_xbasis(i_s - 1);
%     ix_pn_2 = i_xbasis(i_s - 2);
% 
%     Aout(:, i_xbasis(i_s)) = nextPolynomial(Aout(:, ix_pn_1), ...
%         Aout(:, ix_pn_2), n, x);
%     x_order(i_xbasis(i_s)) = n;
%     
%     %y functions
%     iy_pn_1 = i_ybasis(i_s - 1);
%     iy_pn_2 = i_ybasis(i_s - 2);
% 
%     Aout(:, i_ybasis(i_s)) = nextPolynomial(Aout(:, iy_pn_1), ...
%         Aout(:, iy_pn_2), n, y);
%     y_order(i_ybasis(i_s)) = n;
%     
%     mix_q = (i_xbasis(i_s) + 1):(i_ybasis(i_s) - 1);
%     for i_j = 1:numel(mix_q)
%         yo = i_j;
%         xo = n - i_j;
%         xsel = i_xbasis(xo + 1);
%         ysel = i_ybasis(yo + 1);
% 
%         Aout(:, mix_q(i_j)) = Aout(:, xsel) .* Aout(:, ysel);
%         x_order(mix_q(i_j)) = xo;
%         y_order(mix_q(i_j)) = yo;
%     end
% end
% 
% norm_factor = sqrt(2 ./ (2 * x_order + 1)) .* sqrt(2 ./ (2 * y_order + 1));
% for ii = 1:size(Aout,2)
%     Aout(:,ii) = Aout(:,ii) / norm_factor(ii);
% end
% 
% end
% 
% function pn = nextPolynomial(pn_1, pn_2, n, x)
% 
% pn = ((2 * n - 1) * x .* pn_1 - (n - 1) * pn_2) / n;
% 
% end
