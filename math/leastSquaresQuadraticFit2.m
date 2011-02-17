function c = leastSquaresQuadraticFit2(x, y, z)
% c = leastSquaresQuadraticFit2(x, y, z)
% Calculates the least square error best fit quadratic polynomial to data
% collected in R2 -> R.
%
% In other words, this function finds the parameters to fit the polynomial
% z = a x^2 + b xy + c y^2 + dx + ey + f
%
% The parameters a..f are returned in c = [a b c d e f]

if nargin == 1
    z = x(3,:);    
    y = x(2,:);
    x = x(1,:);
end

x = x';
y = y';

A = [x.^2 x.*y y.^2 x y ones(size(x))];

ATA = A' * A;

[U S V] = svd(ATA);

invS = S;
invS(logical(eye(size(S)))) = 1./invS(logical(eye(size(S))));

ATAInv = U * invS * V';

c = (ATAInv * A' * z);