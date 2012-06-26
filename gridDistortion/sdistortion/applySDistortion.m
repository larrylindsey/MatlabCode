function xys = applySDistortion(xy, rk, sk)
% xy - original meshgrid in [n 2]
% rk - radial distortion coefficients, order 1...k
% sk - spiral distortion coefficients, order 0...k-1
%


r = sqrt(sum(xy.^2, 2));
xys = xy;

rdistortion = zeros(size(r));
theta = rdistortion;

for o = 1:numel(rk)
    rdistortion = rdistortion + rk(o) * r.^o;
end

for o = 1:numel(sk)
    theta = sk(o) * r.^(o - 1);
end

for n = 1:size(xy, 1)
    xys(n,:) = xys(n,:) * rotmat(theta(n));
end

xys = xys + xy .* repmat(rdistortion, [1 2]);
