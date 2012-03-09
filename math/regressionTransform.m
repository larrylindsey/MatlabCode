function tr = regressionTransform(ptsin, ptsout, order, type, data)

if nargin < 5
    data.u = [min(ptsin(:,1)) max(ptsin(:,1))];
    data.v = [min(ptsin(:,2)) max(ptsin(:,2))];
    data.n = size(ptsin, 1);
    data.gamma = [];
    data.weight = [];
    if nargin < 4
        type = @legendreMat;
    end
end

tr.order = order;
tr.T = calculateTransMat(ptsin, ptsout, order, data, type);
%tr.Tinv = calculateTransMat(ptsout, ptsin, order, type);
tr.iDim = 2;
tr.oDim = 2;
tr.type = type;
tr.data = data;
tr.Tinv = [];
% tr = populateTransInverse(tr);
