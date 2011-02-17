function [outGrid sliceCollector gridBasis] = ...
    contourGrid(secdoc, posname, negname, r, n, d)


%if nargin > 4
doErode = ~isempty(d) && d > 0;
%else
%    doErode = false;
%end

if ~iscell(posname)
    posname = {posname};
end

if ~iscell(negname) && ~isempty(negname)
    negname = {negname};
elseif isempty(negname)
    negname = {};
end


[gridBasis x y zi dx dy dz] = getGridBasis(n, secdoc, posname);


if doErode
    if ~isinteger(d)
        d = ceil(d / mean([dx dy]));
    end    
    if d <= 0
        doErode = false;
    else
        s = strel('square', d);
    end
end

if isinteger(r)
    r = double(r) * (dx + dy) / 2;
end

sliceCollector = false(size(x, 1), size(x, 2), numel(zi));

outGrid = [];

disp('Done.  Let''s go');

for i_sec = 1:numel(zi)
    sec = zi(i_sec);

    slice = contourSlice(secdoc, posname, sec, x, y);
    
%     slice = false(size(x));
%     for i_p = 1:numel(posname)
%         slice = or(slice, contourSlice(secdoc, posname{i_p}, sec, x, y));
%     end
    
%     for i_n = 1:numel(negname)
%         slice = and(slice, ...
%             not(contourSlice(secdoc, negname{i_n}, sec, x, y)));
%     end

    if doErode
        slice = imerode(slice, s);
    end

    if ~isempty(negname)
        slice = and(slice, not(contourSlice(secdoc, negname, sec, x, y)));
    end    
   
    sliceCollector(:,:,i_sec) = slice;
    
    sliceEdge = mkSliceEdge(slice);
    
    sliceX = x(sliceEdge);
    sliceY = y(sliceEdge);
    
    sliceGrid = quickmesh(cat(2, sliceX, sliceY), r);
    
    sliceX = x(slice);
    sliceY = y(slice);
    
    sliceGrid = quickmesh(cat(2, sliceX, sliceY), r, sliceGrid);
    
    z = repmat(sec * dz, [size(sliceGrid, 1) 1]);
    
    sliceGrid = cat(2, sliceGrid, z);
    
    outGrid = cat(1, outGrid, sliceGrid);
end


end

function sliceEdge = mkSliceEdge(slice)
dSlice1 = cat(1, diff(slice, 1, 1), zeros(1, size(slice, 2)));
dSlice2 = cat(2, diff(slice, 1, 2), zeros(size(slice, 1), 1));
sliceEdge = or(dSlice1 ~= 0, dSlice2 ~= 0);
end

function [gridBasis x y zi dx dy dz] = getGridBasis(n, secdoc, posname)

if isnumeric(n)
    
    if numel(n) == 1
        n = [n n];
    end

    disp('Doing some pre-calculations...');

    xMin = nan;
    xMax = nan;
    yMin = nan;
    yMax = nan;
    ziMin = nan;
    ziMax = nan;

    dz = secdoc(1).section.thickness;

    posContour = cell(1, numel(posname));

    for i_p = 1:numel(posname)
        posContour{i_p} = extractContour(secdoc, posname{i_p});
        if ~isempty(posContour{i_p})
            extPts = cat(1, posContour{i_p}.transPoints);
            sec = [posContour{i_p}.section];
            xMin = min(xMin, min(extPts(:,1)));
            xMax = max(xMax, max(extPts(:,1)));
            yMin = min(yMin, min(extPts(:,2)));
            yMax = max(yMax, max(extPts(:,2)));
            ziMin = min(ziMin, min(sec));
            ziMax = max(ziMax, max(sec));
        end
    end

    dx = (xMax - xMin) / n(1);
    dy = (yMax - yMin) / n(2);


    gridBasis{1} = linspace(xMin, xMax, n(1));
    gridBasis{2} = linspace(yMin, yMax, n(2));

    [x y] = meshgrid(gridBasis{1}, gridBasis{2});
    zi = ziMin:ziMax;

    gridBasis{3} = zi * dz;

    gridBasis{4} = zi;
else
    gridBasis = n;
    [x y] = meshgrid(gridBasis{1}, gridBasis{2});
%    x = gridBasis{1};
%    y = gridBasis{2};
    z = gridBasis{3};
    zi = gridBasis{4};
    
    dx = median(diff(x));
    dy = median(diff(y));
    dz = median(diff(z));
end
end