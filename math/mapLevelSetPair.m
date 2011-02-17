function [map] = mapLevelSetPair(l1, l2)

n1 = size(l1, 1);
n2 = size(l2, 1);

dd = dist2(l1, l2);

ll = cat(1, l1, l2);
vind = 1:size(ll, 1);

vind1 = 1:size(l1, 1);
vind2 = setdiff(vind, vind1);

[junk imin1_2] = min(dd, [], 2);
[junk imin2_1] = min(dd, [], 1);

map = zeros(size(ll, 1), 2);

map(1:numel(vind1), 1) = 1:numel(vind1);
map(1:numel(vind1), 2) = imin1_2;

map(numel(vind1) + (1:numel(vind2)), 1) = imin2_1;
map(numel(vind1) + (1:numel(vind2)), 2) = 1:numel(vind2);

% This should sort map first by the first index in columns, then the
% second, as well as removing redunant mappings.
map = sortmap(map);

%Collect frequency-histograms
freq1 = zeros(size(vind1));
for ii = 1:numel(vind1)
    freq1(ii) = sum(map(:,1) == ii);
end

freq2 = zeros(size(vind2));
for ii = 1:numel(vind2)
    freq2(ii) = sum(map(:,2) == ii);
end

%Find indices that are mapped only once.  Make sure they're part of two
%triangles later on.
singles1 = find(freq1 == 1);
for ii = 1:numel(singles1)
    sloc = find(map(:,1) == singles1(ii));
    pard = map(sloc, 2);
    fwdpard = modi(pard + 1, freq2);
    revpard = modi(pard - 1, freq2);
    pards_locs = find(map(:,2) == pard);
    pards_pards = map(pards_locs,1);
    if freq2(fwdpard) == 1 && ~any(pards_pards > singles1(ii))
        map = cat(1, map, [singles1(ii), fwdpard]);
    elseif freq2(revpard) == 1 && ~any(pards_pards < singles1(ii))
        map = cat(1, map, [singles1(ii), revpard]);
    end
end

map = sortmap(map);

%Recalculate frequency histograms
freq1 = zeros(size(vind1));
for ii = 1:numel(vind1)
    freq1(ii) = sum(map(:,1) == ii);
end

freq2 = zeros(size(vind2));
for ii = 1:numel(vind2)
    freq2(ii) = sum(map(:,2) == ii);
end

singles2 = find(freq2 == 1);
for ii = 1:numel(singles2)
    sloc = find(map(:,2) == singles2(ii));
    pard = map(sloc, 1);
    fwdpard = modi(pard + 1, freq1);
    revpard = modi(pard - 1, freq1);
    pards_locs = find(map(:,1) == pard);
    pards_pards = map(pards_locs,2);
    if freq1(fwdpard) == 1 && ~any(pards_pards > singles2(ii))
        map = cat(1, map, [singles1(ii), fwdpard]);
    elseif freq1(revpard) == 1 && ~any(pards_pards < singles2(ii))
        map = cat(1, map, [singles1(ii), revpard]);
    end
end

map = sortmap(map);


%OK.  Now we have a map that is gauranteed to be nice as long as we have
%convex contours, or reasonably convex contours, but not if one contour is
%much more concave than the other, so we need to fix that.

diffmap1 = mod(diff(map(:,1)), n1);
diffmap2 = mod(diff(map(:,2)), n2);

if any(diffmap1 > 1) || any(diffmap2 > 1)
    theta1 = calculateLSAngle(l1);
    theta2 = calculateLSAngle(l2);

    % First LS
    incpath = calculateIncrementalPaths(map(:,1));
    map = orderByIncrement(map, incpath, theta1, theta2, 1);
    
    % Second LS
    incpath = calculateIncrementalPaths(map(:,2));
    map = orderByIncrement(map, incpath, theta1, theta2, 2);
           
end

%Finally, convert the map to a triangulation
map = sortmap(map);
map = triMap(map);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map = triMap(map)

k = max(map(:,1));
l = max(map(:,2));

n = size(map, 1);
for ii = n:-1:1
    ifwd = modi(ii + 1, n);
    d = map(ifwd, :) - map(ii, :);
    d = [mod(d(1), k) mod(d(2), l)];
    if sum(d) == 2
        newmap = map(ii, :);
        newmap(2) = modi(map(ii,2) + 1, l);
        map = cat(1, map(1:ii, :), newmap, map((ii + 1):end, :));
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta = calculateLSAngle(ls)
n = size(ls, 1);
theta = zeros(n, 1);

for ii = 1:n
    ifwd = modi(ii + 1, n);
    irev = modi(ii - 1, n);
    v = ls(ifwd,:) - ls(irev,:);
    theta(ii) = atan2(v(2), v(1));
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outmap = orderByIncrement(map, incpath, theta1, theta2, index)
if numel(incpath) == 1
    outmap = map;    
elseif index == 2
    outmap = orderByIncrement(map(:,[2 1]), incpath, theta2, theta1, 1);
    outmap = outmap(:, [2 1]);
else
    k = max(map(:,1));
    l = max(map(:,2));
    
    w = zeros(1, numel(incpath));
    for i_p = 1:numel(incpath)
        path = incpath{i_p};
        diffAngle = theta1(map(path, 1)) - theta2(map(path,2));        
        w(i_p) = rms(mod(diffAngle + pi, 2 * pi) - pi);
    end
    [junk iStartPath] = min(w);
    
    outmap = map(incpath{iStartPath}, :);
    
    map = map(setdiff(1:size(map, 1), incpath{iStartPath}), :);
    
    outmap = catMap(outmap, map, k, l);
    
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map = catMap(map, map2, k, l)

mapCropSel = [];

for i_m = 1:size(map2, 1)
    if any(map(:,1) == map2(i_m,1))
        if any(map(:, 2) == map2(i_m, 2))
            mapCropSel = [mapCropSel, i_m];
        end
    end   
end

map2(mapCropSel, :) = [];
noFix = true;

while ~isempty(map2) && noFix
    n = map(end, 1);
    q = mod(map2(:,1) - n, k);
    [d inext] = min(q);
    if d > 1
        noFix = false;
        %error('Assertion Failure');
        disp('wierd');
        
    else

        c = mod(map2(inext, 2) - map(end, 2), l);
        map(end + 1, :) = map2(inext, :);
        if c > 1
            map(end, 2) = map(end - 1, 2);
        end

        map2(inext, :) = [];
    end
end

if ~isempty(map2)
    fixIndices = unique(map2(:, 2));
    for i_m = 1:numel(fixIndices)
        n = fixIndices(i_m);
        q = mod(n - map(:,2), l);
        d = min(q);
        inext = find(q == d, 1, 'last');
        
        newmap = cat(2, map(inext, 1), n);
                
        map = cat(1, map(1:inext,:), newmap, map((inext + 1):end, :));
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function incpath = calculateIncrementalPaths(ind)
%sel = false(numel(ind), 1);
sel = logical([]);
unsel = true(numel(ind), 1);

while any(unsel)
    firstUnsel = find(unsel, 1, 'first');
    ix = ind(firstUnsel);
    unsel(firstUnsel) = false;
    sel = cat(2, sel, false(size(unsel)));
    sel(firstUnsel, end) = true;
    unsel_f = find(unsel);
    for i_f = 1:numel(unsel_f)
        ii = unsel_f(i_f);
        d = ind(ii) - ix;
        if any([d d] == [0 1])
            sel(ii, end) = true;
            unsel(ii) = false;
            ix = ind(ii);
        end
    end
end

incpath = cell(1, size(sel, 2));
for jj = 1:size(sel, 2)
    incpath{jj} = find(sel(:,jj));
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map = sortmap(map)
map = unique(map, 'rows');

n = max(map(:,1));
m = max(map(:,2));

for in = 1:n
    secsel = map(:,1) == in;
    msel = map(secsel,2);
    mseld = diff(msel);
    q = 0;

    while any(mseld > 1) && any(msel == 1) && q < m
        q = q + 1;
        mseld = mod(mseld + q, m) - q;
    end
    
    q = mod(q, m);
    
    if q ~= 0

        [junk isort] = sort(mod(msel + q, m));
        msel = msel(isort);

        map(secsel, 2) = msel;
    end
end

end
