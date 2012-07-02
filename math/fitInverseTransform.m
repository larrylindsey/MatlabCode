function tr = fitInverseTransform(tr, matrixFunction, invParam)
% tr = fitInverseTransform(tr, matrixFunction, invParam)
%   Fits the inverse to the transform tr.
%

if nargin < 3
    invParam.param = tr.param;
    invParam.tpts = [];
    invParam.fpts = [];
    invParam.n = 1024;
    invParam.lim = [-1; 1];
    invParam.order = tr.order;
    
    if nargin < 2
        matrixFunction = tr.matrixFun;
    end
end

fpts = invParam.tpts;
tpts = invParam.fpts;

% If fpts is empty, make our own fpts, tpts.
if isEmpty(fpts)
    lim = invParam.lim;
    n = invParam.n;
    
    if size(lim, 2) == 1
        lim = repmat(lim, [1 d]);
    end
    if numel(n) == 1
        n = repmat(round(n^(1/d)), [d 1]);
    end

    t = cell(1,d);
    
    for ii = 1:d
        t{ii} = linspace(lim(1,ii), lim(2,ii), n(ii));
    end
    
    fpts = gridRC(t{:});
    
    tpts = applyTransform(fpts, tr);    
end

tr.inv = fitTransform(fpts, tpts, invParam.order, matrixFunction, ...
    invParam.param);

end
