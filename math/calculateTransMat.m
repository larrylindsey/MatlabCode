function tr = calculateTransMat(fpts, tpts, order, type)

if nargin < 4
    type = @legendreMat;
end

A = type(fpts(:,1), fpts(:,2), order);

tr = cat(2, leastSquaresFit(A, tpts(:,1)), leastSquaresFit(A, tpts(:,2)));

if size(tr, 1) <= 1
    tr = cat(1, tr, eye(2));
end