function eaOutput = analyzeSiftCacheAlignment(siftCache, tr)
% eaOutput = analyzeSiftCacheAlignment(siftCache, tr, transFun)
%  Measures the alignment error of a 
%
%


if nargin < 2
    tr = identityTransform();
end


if ischar(siftCache)
    siftCache = load(siftCache);
end

n = size(siftCache.feat1, 2);

eaOutput.sectionWise = cell(n, 1);
eaOutput.d = cell(n, 1);

for ii = 1:n
    xy1 = siftCache.feat1{ii};
    xy2 = siftCache.feat2{ii};
    
    xy1tr = applyTransform(xy1, tr);
    xy2tr = applyTransform(xy2, tr);
    
    eaOutput.d{ii} = sqrt(sum((xy1).^2,2));
    
    if isempty(xy1)
        eaOutput.sectionWise{ii} = [];
    else
        eaOutput.sectionWise{ii} = measure(xy1tr, xy2tr);
    end
end

eaOutput.all = cat(1, eaOutput.sectionWise{:});

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = measure(xy1, xy2)

xy2aff = affineAlign(xy2, xy1);

e = sqrt(sum((xy2aff - xy1).^2,2));
end




