function c = dimWiseLeastSquaresFit(A, z, param)
% c = dimWiseLeastSquaresFit(A, z, <param>)
%   
%

if nargin < 1
    c = leastSquaresFit;
    return;
end

d = size(z, 2);
c = leastSquaresFit(A(:,:,1), z(:,1), param);
c = repmat(c(:), [1 d]);

for ii = 2:d
    c(:,d) = leastSquaresFit(A(:,:,d), z(:,d), param);
end

