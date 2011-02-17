function [tri tripts]  = triangulateContour(secdoc, cname)

if nargin > 1
    contour = extractContour(secdoc, cname);
else
    contour = secdoc;
end

for i_c = 1:numel(contour)
    contour(i_c).centroid = contourCentroid(contour(i_c));
end

zall = [contour.z];
zuid = sort(unique(zall));

tri = [];
tripts = [];

for i_z = 2:numel(zuid)
    z1 = zuid(i_z - 1);
    z2 = zuid(i_z);
    
    slice1 = contour(zall == z1);
    slice2 = contour(zall == z2);
    
    [pts1 pts2] = mapPoints(slice1, slice2);
    
    for i_p = 1:numel(pts1)
        [t p] = levelSetMesh({pts1{i_p}, pts2{i_p}}, [z1 z2]);
        [tri tripts] = catTri(tri, tripts, t, p);
    end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yn = intersection(slice1, slice2)
pts1 = slice1.transPoints;
pts2 = slice2.transPoints;

minrange = max(min(pts1, [], 1), min(pts2, [], 1));
maxrange = min(max(pts1, [], 1), max(pts2, [], 1));

if any(minrange > maxrange)
    yn = false;
    return;
end

[x y] = meshgrid(linspace(minrange(1), maxrange(1), 64), ...
    linspace(minrange(2), maxrange(2), 64));

boolslice1 = inpolygon(x, y, pts1(:,1), pts1(:,2));
boolslice2 = inpolygon(x, y, pts2(:,1), pts2(:,2));

yn = any(boolslice1 & boolslice2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pts = doMergePts(contours)
if numel(contours) == 1
    pts = contours.transPoints;
else
    pts = mergeContour(contours(1).transPoints,...
        contours(2).transPoints);
    for i_c = 3:numel(contours)
        pts = mergeContour(pts, contours(i_c).transPoints);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pts1 pts2] = mapPoints(slice1, slice2)

ctr1 = cat(1, slice1.centroid);
ctr2 = cat(1, slice2.centroid);
 
dmap = dist2(ctr1, ctr2);

 
[junk imin1] = min(dmap, [], 2);
[junk imin2] = min(dmap, [], 1);


sliceCat = [slice1 slice2];
idx = [imin1(:)'+numel(imin1) imin2(:)'];
idf = zeros(size(idx)) + 2;
idf(1:numel(slice1)) = 1;


for ii = 1:numel(idx)
    if ~intersection(sliceCat(ii), sliceCat(idx(ii)))
        idx(ii) = ii;
    end
end


cycles = groupByCycle(idx);
ct = 1;
for i_c = 1:numel(cycles)
    cycle = cycles{i_c};
    if numel(cycle) > 1
        i1 = cycle(idf(cycle) == 1);
        i2 = cycle(idf(cycle) == 2);
        
        pts1{ct} = doMergePts(sliceCat(i1));
        pts2{ct} = doMergePts(sliceCat(i2));
        ct = ct + 1;
    end    
end

% pts1 = {};
% pts2 = {};
% 
% for ii = 1:numel(slice1)
%     pts1 = {pts1{:} slice1(ii).transPoints};
% end
% 
% for ii = 1:numel(slice2)
%     pts2 = {pts2{:} slice2(ii).transPoints};
% end
% 
% imin1u = unique(imin1);
% for ii = 1:numel(imin1u)
%     mergeIndices = find(imin1 == imin1u(ii));
%     if numel(mergeIndices) > 1
%         pushIndex = mergeIndices(1);
%         pullIndices = mergeIndices(2:end);
%         mergePts = pts1{pushIndex};
%         for imerge = 1:numel(pullIndices)
%             mergePts = mergeContourMesh(mergePts,...
%                 pts1{pullIndices(imerge)});
%             imin2(imin2 == pullIndices(imerge)) = pushIndex;
%         end
%         pts1{pushIndex} = mergePts;
%     end
% end
% 
% imin2u = unique(imin2);
% for ii = 1:numel(imin2u)
%     mergeIndices = find(imin2 == imin2u(ii));
%     if numel(mergeIndices) > 1
%         pushIndex = mergeIndices(1);
%         pullIndices = mergeIndices(2:end);
%         mergePts = pts2{pushIndex};
%         for imerge = 1:numel(pullIndices)
%             mergePts = mergeContourMesh(mergePts,...
%                 pts2{pullIndices(imerge)});
%             imin1(imin1 == pullIndices(imerge)) = pushIndex;
%         end
%         pts2{pushIndex} = mergePts;
%     end
% end
% 
% pts1 = {pts1{sort(unique(imin2))}};
% pts2 = {pts2{sort(unique(imin1))}};

%keyboard
end
