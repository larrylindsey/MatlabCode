function trblock = makeTrBlock(fpts, tpts, paramblock)

trblock = fitTransform(fpts, tpts, paramblock(1).order, paramblock(1).fun);
trblock = repmat(trblock, size(paramblock));

for ii = 2:numel(paramblock)
    trblock(ii) = fitTransform(fpts, tpts, paramblock(ii).order, ...
        paramblock(ii).fun);
end

trblock = setDefaultData(trblock);

end
