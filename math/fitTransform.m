function tr = fitTransform(fpts, tpts, order, matrixFunction, param)
% tr = fitTransform(fpts, tpts, order, [matrixFunction, [param]])
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
% matrixFunction [function handle] - the handle to a function that takes
%        arguments (tpts, fpts, order), and returns a single output, A.
%        This is typically a matrix used in a least squares calculation.
%        When left out, defaults to @legendreMat
% param [1] struct - an optional parameter struct that is passed to the fit
%        function. If left out, it is assigned the output of
%        createTransformFun(), with no arguments.

if nargin < 4
    matrixFunction = @legendreMat;
end

[doTransformFun, createTransformFun] = matrixFunction();

if nargin < 5
    param = createTransformFun();
end

if ~isfield(param, 'fun')
    param.fun = @taylorMat;
end

[A mOrder] = matrixFunction(fpts, tpts, order, param);

T = createTransformFun(A, tpts, param);

tr.T = T;
tr.inv = [];
tr.doTrans = doTransformFun;
tr.createTrans = createTransformFun;
tr.matrixFun = matrixFunction;
tr.param = param;
tr.order = order;
tr.monomialOrder = mOrder;
tr.ndim = size(fpts, 2);

end
