function eaOutput = analyzeSiftCacheAlignment(siftCache, tr)

if ischar(siftCache)
    siftCache = load(siftCache);
end

n = size(siftCache.feat1, 2);

eaOutput.raw = cell(n, 1);
eaOutput.trans = cell(n, 1);

for ii = 1:n
    xy1 = siftCache.feat1{ii};
    xy2 = siftCache.feat2{ii};
    
    xy1tr = doTransform(xy1, tr);
    xy2tr = doTransform(xy2, tr);
    
    if isempty(xy1)
        eaOutput.raw{ii} = inf;
        eaOutput.trans{ii} = inf;
    else
        eaOutput.raw{ii} = measure(xy1, xy2);
        eaOutput.trans{ii} = measure(xy1tr, xy2tr);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = measure(xy1, xy2)

xy2aff = affineAlign(xy2, xy1);

e = sqrt(mean((xy2aff - xy1).^2,2));
end




