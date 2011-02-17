function [grid ec] = conditionMergeGrids(sliceblock, gridbasis, varargin)

if mod(numel(varargin), 2) == 1
    extContour = varargin{1};
    varargin = {varargin{2:end}};
else
    extContour = [];
end

if numel(varargin) == 2
    grid = varargin{1};       
    w = varargin{2};
    if ~iscell(grid)
        grid = {grid};
    end
else
    grid = {varargin{1:2:end}};
    w = [varargin{2:2:end}];
end

ng = zeros(numel(grid), 1);
gb = ng;
ge = ng;

for i_g = 1:numel(grid)
    ng(i_g) = size(grid{i_g}, 1);
    if i_g == 1
        gb(i_g) = 1;
    else
        gb(i_g) = ge(i_g - 1) + 1;
    end
    ge(i_g) = gb(i_g) + ng(i_g) - 1;
end

grid = cat(1, grid{:});

de = delaunay3(grid(:,1), grid(:,2), grid(:,3));

I = de(:);
J = de(:,[2 3 4 1]); J = J(:);

[I J] = conditionEdges(sliceblock, I, J, grid, gridbasis, extContour);

grpI = groupEdges(I, gb, ge);
grpJ = groupEdges(J, gb, ge);

w = max(cat(1, w(grpI), w(grpJ)), [], 1)';

D = rms(grid(I,:) - grid(J,:), 2) .* w;

ec = cat(2, I, J, D);

end

function [I J] = conditionEdges(sliceblock, I, J, grid, gridbasis, ...
    extContour)

dx = median(diff(gridbasis{1}));
dy = median(diff(gridbasis{2}));
dz = median(diff(gridbasis{3}));

x0 = gridbasis{1}(1);
y0 = gridbasis{2}(1);
z0 = gridbasis{3}(1);

spt = grid(I,:);
fpt = grid(J,:);

six = (spt(:,1) - x0) / dx;
siy = (spt(:,2) - y0) / dy;
siz = (spt(:,3) - z0) / dz;

fix = (fpt(:,1) - x0) / dx;
fiy = (fpt(:,2) - y0) / dy;
fiz = (fpt(:,3) - z0) / dz;

zReject = logical(abs(fiz - siz) >= 1.5);

D = rms(fpt - spt, 2);


sameSliceSel = logical(abs(fiz - siz) < .5);

if isempty(extContour)
    cx = round((six + fix) * .5);
    cy = round((siy + fiy) * .5);
    cz = round((siz + fiz) * .5);

    cc = linearIndex(size(sliceblock), cat(2, cx, cy, cz));

    ccok = (cc > 0) & ( cc <= numel(sliceblock));

    cctest = true(size(cc));

    cctest(ccok) = sliceblock(cc(ccok));

    xyReject = and(sameSliceSel, cctest);


    dReject = D > 2.5 * median(D);

    reject = xyReject | zReject | dReject;
else
    czdiff = diff(sort([extContour.z]));
    czdiff = czdiff(czdiff ~= 0);
    dz = median(czdiff);
    
    dTest = D > sqrt(2) * min(D);
    
    sameSliceSel = sameSliceSel(dTest);
    sameSliceSel = find(sameSliceSel);
    xyReject = false(size(sameSliceSel));
    z = (spt(sameSliceSel,3) + fpt(sameSliceSel,3)) / 2;
    for i_slice = 1:numel(sameSliceSel)
        index = sameSliceSel(i_slice);
        gridspt = grid(I(index),1:2);
        gridfpt = grid(J(index),1:2);
        l1 = cat(1, gridspt, gridfpt)';
        
        testContourIndices = min(abs([extContour.z] - z(i_slice)))...
            < (dz / 2);
        
        testContour = extContour(testContourIndices);
                      
        
        for ilset = 1:numel(testContour)

            lset = testContour(ilset).transPoints;

            i_pt = 1;
            while i_pt < size(lset, 1) && ~xyReject(i_slice)
                cspt = lset(i_pt,:);
                cfpt = lset(i_pt + 1, :);
                l2 = cat(1, cspt, cfpt)';
                [junk xyReject(i_slice)] = intersection(l1, l2);
                i_pt = i_pt + 1;
            end

            if ~xyReject(i_slice)
                cspt = lset(end,:);
                cfpt = lset(1, :);
                l2 = cat(1, cspt, cfpt)';
                [junk xyReject(i_slice)] = intersection(l1, l2);
            end
        end
        
    end
    
    reject = zReject;
    reject(sameSliceSel) = reject(sameSliceSel) | xyReject;
end

I(reject) = [];
J(reject) = [];

end

function grp = groupEdges(I, gb, ge)
grp = zeros(size(I));
for i_g = 1:numel(gb)
    sel = logical(I >= gb(i_g) & I <= ge(i_g));
    grp(sel) = i_g;
end
end
