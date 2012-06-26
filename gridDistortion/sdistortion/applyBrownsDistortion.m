function xys = applyBrownsDistortion(xy, k, p)
% xy - original meshgrid in [n 2]
% rk - radial distortion coefficients, order 1...k
% sk - spiral distortion coefficients, order 0...k-1
%


r2 = sum(xy.^2, 2);

r_rep = repmat(sqrt(r2), [1 2]);


krpoly = zeros(size(r_rep));
for ii = 1:numel(k)
    krpoly = krpoly + k(ii)*(r_rep.^ii);
end

prpolyx = p(1) * (r2 + 2*(xy(:,1).^2)) + 2 * p(2) * (xy(:,1).*xy(:,2));
prpolyx = prpolyx .* (1 + p(3)*r2);

prpolyy = p(2) * (r2 + 2*(xy(:,1).^2)) + 2 * p(1) * (xy(:,1).*xy(:,2));
prpolyy = prpolyy .* (1 + p(3)*r2);

prpoly = cat(2, prpolyx, prpolyy);

xys = xy + xy.*krpoly + prpoly;
