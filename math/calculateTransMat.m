function tr = calculateTransMat(fpts, tpts, order, param, type)

if nargin < 5
    type = @legendreMat;
    if nargin < 4
        param.gamma = [];
        param.weight = [];
    end
end

A = type(fpts(:,1), fpts(:,2), order);

tr = cat(2, leastSquaresFit(A, tpts(:,1), param), ...
    leastSquaresFit(A, tpts(:,2), param));

if size(tr, 1) <= 1
    tr = cat(1, tr, eye(2));
end