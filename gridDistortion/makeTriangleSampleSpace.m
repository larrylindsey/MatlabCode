function [sampleSpace angles lvectexpect inds] = makeTriangleSampleSpace(rc,...
    expectAngle)

if nargin < 2
    expectAngle = pi / 2;
end


%setup basic variables
rotmat = [cos(expectAngle) -sin(expectAngle);...
    sin(expectAngle) cos(expectAngle)];
s = size(rc, 1);
iid = 1:s;

%find pairwise closest vectors
dctr = dist2(rc, rc);
%set auto-distance to infinity, otherwise we find that each vector is
%closest to itself.
dctr(sub2ind([s s], 1:s, 1:s)) = inf;
[junk imind] = min(dctr);

%rvect contains the "right hand" vector
rvect = rc(imind, :);
d_rvect = rvect - rc;

%lvectexpect is the location of the expected "left hand" vector
lvectexpect = (rotmat * d_rvect')' + rc;
dexpect = dist2(lvectexpect, rc);
dexpect(sub2ind([s s], 1:s, 1:s)) = inf;
%find the vector that is closest to the expected "left hand" location.
[junk imine] = min(dexpect, [], 2);

% Indicates indices into rc, for each triangle found.
inds = cat(2, imind(:), iid(:), imine(:));

% Find the actual "left hand" vector
lvect = rc(imine, :);
d_lvect = lvect - rc;

% Normalize the offset vectors.
d_lvect_norm = sqrt(dot(d_lvect, d_lvect, 2));
d_rvect_norm = sqrt(dot(d_rvect, d_rvect, 2));
nd_lvect = d_lvect ./ repmat(d_lvect_norm, [1 2]);
nd_rvect = d_rvect ./ repmat(d_rvect_norm, [1 2]);

% Calculate angles
d_cos = dot(nd_rvect, nd_lvect, 2);
angles = acos(d_cos);

%%%% Collate / Orient Vectors %%%%
% Generally, vectors declared right- or left-handed should point in the
% same direction.

% Reference pair should be the one that is closest to 90-degrees.
[junk isqr] = min(abs(d_cos));

lref = nd_lvect(isqr,:);
rref = nd_rvect(isqr,:);

%Check left vect dot product
lldot = nd_lvect * lref';
lrdot = nd_lvect * rref';

% If a "left-handed" offset vector matches the right-hand reference better
% than the left-hand reference, it is really a right-handed vector, and
% should be swapped with its sibling.
flipsel = abs(lrdot) > abs(lldot);

d_temp = d_lvect;
d_lvect(flipsel,:) = d_rvect(flipsel,:);
d_rvect(flipsel,:) = d_temp(flipsel,:);

% Recalculate rather than do tedious bookkeeping
d_lvect_norm = sqrt(dot(d_lvect, d_lvect, 2));
d_rvect_norm = sqrt(dot(d_rvect, d_rvect, 2));
nd_lvect = d_lvect ./ repmat(d_lvect_norm, [1 2]);
nd_rvect = d_rvect ./ repmat(d_rvect_norm, [1 2]);

% Now the ith vector in d_lvect is approximately orthogonal to the ith
% vector in d_rvect, however, many of these vectors will point in the
% opposite direction of their nominal reference.  Find those vectors, and
% negate them:
lldot = nd_lvect * lref';
rrdot = nd_rvect * rref';

l_neg_sel = lldot < 0;
r_neg_sel = rrdot < 0;

d_lvect(l_neg_sel,:) = -d_lvect(l_neg_sel,:);
d_rvect(r_neg_sel,:) = -d_rvect(r_neg_sel,:);

% Finally, for each location in rc, find the indices in rc for those
% vectors that are within a neighborhood of it.

nbd_dist = mean([d_lvect_norm(isqr) d_rvect_norm(isqr)]) * 2;
[nbd_r nbd_c] = find(dctr <= nbd_dist^2);
maxsz_nbd = 0;
for ii = 1:size(rc, 1)
    sel = nbd_r == ii;
    nbd{ii} = nbd_c(sel);%#ok
    if numel(nbd{ii}) > maxsz_nbd
        maxsz_nbd = numel(nbd{ii});
    end
end
nbd_block = zeros(size(rc, 1), maxsz_nbd) - 1;
for ii = 1:size(rc, 1)
    locnbd = nbd{ii};
    nbd_block(ii,1:numel(locnbd)) = locnbd;
end

sampleSpace = cat(2, d_lvect, d_rvect, iid(:), nbd_block);
