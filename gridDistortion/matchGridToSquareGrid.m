function [rc_grid] = matchGridToSquareGrid(rc, grid_model, t)
ngrp = 4;

n = size(rc, 1);

% cvect in [n 5]
% cvect(:,1:2) - edge connection indices for rc. IE, let k : 1 <= k <= m,
%                where m is the number of rows is cvect.
%                Let i1 = cvect(k,1), i2 = cvect(k,2), then there exists an
%                edge between the points rc(i1,:) and rc(i2,:).
% cvect(:,3) - for each edge as above, indicates the error associated with
%              it.
% cvect(:,4:5) - grid_model coordinates corresponding to the difference in
%                locations between vertices defining the edge.
%  
cvect = makeLinks(rc, grid_model, [1 0; 0 1; 2 0; 0 2; 1 1; 2 2; 3 0; 0 3]');

% Select edges with error less than t
sel = cvect(:,3) < t;

% Remove redundant edges. Might as well keep the order, too.
cvect = cvect(sel,:);
[junk iu] = unique(cvect(:,1:2), 'rows');%#ok
cvect = cvect(sort(iu),:);

% Separate cvect into an edge list (conn), and a grid model coordinate list
% (conncoord).
conn = cvect(:,1:2);
% cdist = cvect(:,3);
conncoord = cvect(:,4:5);

% Let 1 <= j,k <= n
% The locations rc(j,:), and rc(k,:) define an edge if lm(j,k) == 1
% In this case, rc(k,:) is located at a position approximated by
% rc(j,:) + grid_model(:,1) * rConnMat(j,k) + 
%       grid_model(:,2) * cConnMat(j,K)
[linkMat rConnMat cConnMat] = makeSparseMats(conn, conncoord, n);

% ci contains an index for each row in rc, indicating the component index
% of that point. ciOrder contains a list of the component indices in
% descending order of cardinality.
[ci ciOrder cis] = getComponents(linkMat);

ngrp = min(ngrp, numel(ciOrder));
joingrp = find((8 * cis) > cis(1), 1, 'last');

grp = cell(1, ngrp);

% Mesh the groups, ie, assign each one a row and column
for ig = 1:ngrp
    grp{ig} = meshGroup(ci, ciOrder(ig), rConnMat, cConnMat, linkMat);
end

bigGroup = grp{1};

for ig = 2:joingrp
    bigGroup = joinGroups(bigGroup, grp{ig}, rc, grid_model);
end
    
%[rc_grid c] = groupDist(rc, ci, icisort(1), icisort(2), grid_model);
rc_grid = bigGroup;
rc_grid = rc_grid(:,[1 3 2]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lm rcm ccm] = makeSparseMats(conn, conncoord, n)
lm = sparse([conn(:,1) conn(:,2)], [conn(:,2) conn(:,1)], ...
    ones(2 * size(conn,1), 1) , n, n);
rcm = sparse([conn(:,1) conn(:,2)], [conn(:,2) conn(:,1)], ...
    cat(1, -conncoord(:,1), conncoord(:,1)) , n, n);
ccm = sparse([conn(:,1) conn(:,2)], [conn(:,2) conn(:,1)], ...
    cat(1, -conncoord(:,2), conncoord(:,2)) , n, n);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grpCat = joinGroups(grp1, grp2, rc, gm)

grc1 = rc(grp1(:,1),:);
grc2 = rc(grp2(:,1),:);

dd = dist2(grc1, grc2);
ddminr = repmat(min(dd, [], 1), [size(dd, 1), 1]);
ddminc = repmat(min(dd, [], 2), [1, size(dd, 2)]);
minmin = and((ddminr == dd) , (ddminc == dd));

[minsel1 minsel2] = find(minmin);

gminrc1 = grc1(minsel1,:);
gminrc2 = grc2(minsel2,:);

gmincoord1 = grp1(minsel1,2:3);
gmincoord2 = grp2(minsel2,2:3);

diffrc = gminrc2 - gminrc1;

dcoord = round(diffrc / (gm'));

bdiff = gmincoord1 + dcoord - gmincoord2;

bdiffMedian = median(bdiff, 1);

errno = sum(dist2(bdiff, bdiffMedian) > 0);
fprintf('Error fraction: %d of %d\n', errno, size(bdiff,1));

grp2(:,2) = grp2(:,2) + bdiffMedian(1);
grp2(:,3) = grp2(:,3) + bdiffMedian(2);

grpCat = cat(1, grp1, grp2);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grp = meshGroup(ci, ici, rcm, ccm, lm)
grpList = find(ci == ici);
linked = false(size(grpList));

grpLM = lm(grpList,grpList);

coord = zeros(numel(grpList), 2);
linked(1) = true;
q = 1;

while ~isempty(q)
    currNode = q(1);
    linkNodes = find(and(grpLM(:,currNode), not(linked)));
    
    
    if numel(q) == 1
        q = linkNodes(:);
    else
        q = cat(1, q(2:end), linkNodes(:));
    end
    
    if ~isempty(q)       
        linkSel = grpList(linkNodes);
        coord(linkNodes,1) = rcm(linkSel, grpList(currNode)) + ...
            coord(currNode, 1);
        coord(linkNodes,2) = ccm(linkSel, grpList(currNode)) + ...
            coord(currNode, 2);
        
        linked(linkNodes) = true;
    end
end

if any(not(linked))
    error('That wasn''t supposed to happen');
end

grp = cat(2, grpList, coord);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ci ciOrder cis] = getComponents(lm)
[ci cis] = components(lm);
[cis ciOrder] = sort(cis, 'descend');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cx = makeLinks(rc, grid_model, coord)

cx = [];
for icoord = 1:size(coord, 2)
  cvect = getLinks(rc, grid_model, coord(:,icoord));
  cvect = cat(2, cvect, repmat(coord(:,icoord)', [size(cvect, 1), 1]));
  cx = cat(1, cx, cvect);
end

sw = cx(:,1) > cx(:,2);
cx(sw,1:2) = cx(sw,[2 1]);
cx(sw,4:5) = -cx(sw, 4:5);

% cvect = getLinks(rc, grid_model(1,:), 1);
% cvect = cat(1, cvect, getLinks(rc, grid_model(2,:), 1));
% cvect = cat(1, cvect, getLinks(rc, grid_model(1,:), 2));
% cvect = cat(1, cvect, getLinks(rc, grid_model(2,:), 2));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cvect = getLinks(rc, gm, f)

n = size(rc, 1);

e2 = abs(det(gm));

cvect = zeros(n, 3);
ok = false(n,1);
% sw = false(n,1);

ovect = f' * gm;

for irc = 1:n
    vExpect = rc(irc, :) + ovect;
    dists = dist2(vExpect, rc);
    %[dists iv] = sort(dists, 'ascend');
    [distMin iv] = min(dists);
    
%     ok(irc) = dists(1) < e2 && iv(1) ~= irc;
%     sw(irc) = irc > iv(1);
%     cvect(irc, :) = [irc iv(1) sqrt(dists(1))];

    ok(irc) = distMin < e2 && iv ~= irc;
    %sw(irc) = irc > iv;
    cvect(irc, :) = [irc iv sqrt(distMin)];


end

%cvect(sw,:) = cvect(sw, [2 1 3]);
cvect = cvect(ok,:);

end
