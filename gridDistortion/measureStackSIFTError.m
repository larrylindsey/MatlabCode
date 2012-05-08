function e = measureStackSIFTError(trArray)
e.delta = [];
e.rms = [];
e = repmat(e, [1 numel(trArray)]);
for ii = 1:numel(trArray)
    tr = trArray{ii};
    if isfield(tr, 'T') && ~isempty(tr.T)
        fpts = tr.fromPts;
        tpts = tr.toPts;
        fptsTR = doTransform(fpts, (tr));
        delta = sqrt(sum((fptsTR - tpts).^2, 2));
        e(ii).delta = delta;
        e(ii).rms = rms(delta);
    end
    
end
end