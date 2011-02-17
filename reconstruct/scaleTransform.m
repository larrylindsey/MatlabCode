function tdata = scaleTransform(tdata, factor)

tdataext = tdata(11:end,:);
tdata = tdata(1:10,:);

if ~isempty(tdataext)
    tdataext = scaleTransform(tdataext, factor);
end

scalemat = cat(1, repmat(factor^2, [4 1]), repmat(factor, [3 1]),...
    repmat(1, [2 1]), 1 / factor);
scalemat = repmat(scalemat, [1 2]);
tdata = tdata .* scalemat;

tdata = cat(1, tdata, tdataext);

end