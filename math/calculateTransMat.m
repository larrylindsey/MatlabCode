function tr = calculateTransMat(fpts, tpts, order, param, matrixFunction, ...
    createTransformFun)
% tr = calculateTransMat(fpts, tpts, order, [param, [matrixFunction,
%      [createTransformFun]]])
%   Creates a transform struct by fitting fpts to tpts
%
% fpts [k n] - "from points," k locations in n dimensions.
% tpts [k n] - "to points," k locations in n dimensions.
%        For each i, 0 < i <= k, fpts(i,:) must correspond to
%        tpts(i,:). tr will represents a transform T(x) = y such
%        that T(fpts(i,:) approximately equals T(tpts(i,:))
%        Note: in practice, n usually is 2.
% order [1] - the order of the transform. Most of the time, this represents
%        the maximal polynomial order, and tr is a polynomial transform.
% param [1] struct - an optional parameter struct that is passed to the fit
%        function. If left out, it is assigned the output of
%        createTransformFun(), with no arguments.
% matrixFunction [function handle] - the handle to a function that takes
%        arguments (tpts, fpts, order), and returns a single output, A.
%        This is typically a matrix used in a least squares calculation.
%        When left out, defaults to @legendreMat
% createTransformFun [function handle] - the handle to a function that
%        takes arguments (A, tpts, fpts, param), where A is the output from
%        matrixFunction.
%        When left out, defaults to @leastSquaresFit.
%   


if nargin < 6
    createTransformFun = @leastSquaresFit;
    if nargin < 5
        matrixFunction = @legendreMat;
        if nargin < 4
            param = createTransformFun();
        end
    end
end

A = matrixFunction(fpts, tpts, order);

tr = createTransformFun(A, tpts, fpts, param);

if size(tr, 1) <= 1
    tr = cat(1, tr, eye(2));
end
end
