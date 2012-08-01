function tr = seriesAlignmentTransforms(cache, idIndex, preTr)

doPreTr = nargin > 2 && ~isempty(preTr);
if ~doPreTr
    preTr = [];
end
   
tr = identityTransform(@taylorMat, 1);
n = numel(cache.f);

tr = repmat(tr, [n 1]);

for ii = (idIndex + 1):n
    ptsTo = cache.feat1{ii - 1};
    ptsFrom = cache.feat2{ii - 1};
    
    if doPreTr
        ptsTo = applyTransform(ptsTo, preTr);
        ptsFrom = applyTransform(ptsFrom, preTr);
    end
    
    if isempty(ptsTo)
        tr(ii) = tr(ii - 1);
    else        
        ptsTo = applyTransform(ptsTo, tr(ii - 1));        
        tr(ii) = fitTransform(ptsFrom, ptsTo, 1, @taylorMat);        
    end
    
end


for ii = (idIndex - 1):-1:1
    ptsTo = cache.feat2{ii};
    ptsFrom = cache.feat1{ii};
    
    if doPreTr
        ptsTo = applyTransform(ptsTo, preTr);
        ptsFrom = applyTransform(ptsFrom, preTr);
    end
    
    if isempty(ptsTo)
        tr(ii) = tr(ii + 1);
    else
        ptsTo = applyTransform(ptsTo, tr(ii + 1));
        tr(ii) = fitTransform(ptsFrom, ptsTo, 1, @taylorMat);
    end
end


for ii = 1:numel(tr)
    tr(ii) = fitInverseTransform(tr(ii));
    tr(ii).preTr = preTr;
end

tr = setDefaultData(tr);
end
