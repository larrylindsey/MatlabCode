function [tr outsel] = ransacRegressionTransform(rparams, ptsin, ptsout, order,...
    type, data)

if nargin < 1
    rparams.n = 16;
    rparams.metric = [];
    rparams.maxError = .001;
    rparams.minInliers = 64;
    rparams.maxIter = 1000;
    tr = rparams;
    return;
end

if ~all(size(ptsin) == size(ptsout))
    error('ptsin must have the same size as ptsout');
end


if nargin < 6
    data.u = [min(ptsin(:,1)) max(ptsin(:,1))];
    data.v = [min(ptsin(:,2)) max(ptsin(:,2))];
    data.n = size(ptsin, 1);
    if nargin < 5
        type = @legendreMat;
    end
end

if isempty(rparams.metric)
    rparams.metric = @errorMetric;
end
    

extra.ptsin = ptsin;
extra.ptsout = ptsout;
extra.order = order;
extra.type = type;
extra.data = data;

universe = (1:size(ptsin, 1))';

[outsel, tr] = ransac(universe, @makeModel, rparams.metric, rparams.maxError,...
    rparams.minInliers, rparams.n, rparams.maxIter, extra);
end

function tr = makeModel(sel, extra)
tr = regressionTransform(extra.ptsin(sel,:), extra.ptsout(sel,:), ...
    extra.order, extra.type, extra.data);
end

function measure = errorMetric(tr, sel, extra)
ptsin = extra.ptsin(sel,:);
ptsout = extra.ptsout(sel,:);
ptstr = doTransform(ptsin, tr);

measure = rms(sqrt(sum((ptsout - ptstr).^2, 2)));
end

