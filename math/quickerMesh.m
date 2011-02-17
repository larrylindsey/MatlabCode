function pts = quickerMesh(inPts, r, seedpts)

inPts = cat(1, seedPts, inPts);

sel = true([size(inPts, 1), 1]);
ii = 1;

while ~isempty(ii) && ii <= size(inPts, 1)
    dd = dist2(inPts(ii,:), inPts);
    crud = logical(dd <= r);
    crud(ii) = false;
    sel(crud) = false;
    selh = sel;
    selh(1:ii) = false;
    ii = find(selh, 1, 'first');
end

sel(1:size(seedpts, 1)) = true;

pts = inPts(sel,:);