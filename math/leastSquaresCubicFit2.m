function [c A]= leastSquaresCubicFit2(x, y, z)
% c = leastSquaresCubicFit2(x, y, z)
% Calculates the least square error best fit cubic polynomial to data
% collected in R2 -> R.
%
% In other words, this function finds the parameters to fit the polynomial
% z = a x^3 + b x^2 * y + c x y^2 + d y^3 + 
%     e x^2 + f x * y + g y^2 + 
%     h x + i y +
%     j
%
% The parameters a..j are returned in c = [a b c d e f g h i j]


if nargin == 1
    z = x(3,:);    
    y = x(2,:);
    x = x(1,:);
end

if size(x, 1) < size(x, 2)
    x = x';
end

if size(y, 1) < size(y, 2)
    y = y';
end

if size(z, 1) < size(z, 2)
    z = z';
end

x2 = x.^2;
y2 = y.^2;
x3 = x.^3;
y3 = y.^3;

A = [x3 x2.*y x.*y2 y3 x2 x.*y y2 x y ones(size(x))];

ATA = A' * A;

[U S V] = svd(ATA);

invS = S;
invS(logical(eye(size(S)))) = 1./invS(logical(eye(size(S))));

ATAInv = U * invS * V';

c = (ATAInv * A' * z);