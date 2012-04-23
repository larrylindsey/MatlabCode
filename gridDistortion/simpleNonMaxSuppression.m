function rc = simpleNonMaxSuppression(rc, crossEnergy, d)

d = d.*d;

rcEnergy = crossEnergy(sub2ind(size(crossEnergy), rc(:,1), rc(:,2)));
sel = true(size(rcEnergy));

[~, iSort] = sort(rcEnergy, 'descend');
rc = rc(iSort,:);

currI = 1;

while (numel(currI) > 0 && currI < numel(sel))    
    checkSel = sel;
    checkSel(1:currI) = false;
    checkSelI = find(checkSel);            
    d2 = dist2(rc(currI,:), rc(checkSel,:));
    sel(checkSelI(d2 < d)) = false;
    currI = find(and(sel, checkSel), 1, 'first');
end

rc = rc(sel,:);

end
