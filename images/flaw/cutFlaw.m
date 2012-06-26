function [outPath, cacheStr] = cutFlaw(im, cacheStr)
%

%%% params %%%
if nargin < 2
    cacheStr.n_varblock = 16;
    cacheStr.o_varblock = 0.5;
    cacheStr.n_pathblock = 16;
    cacheStr.b_pathblock = 10;   
end

[n_varblock, o_varblock, n_pathblock, b_pathblock] = parseParms(cacheStr);

if nargin < 1
    outPath = cacheStr;
    return;
end

im = im2single(im(:,:,1));


% Build the flaw map tile image
if ~isfield(cacheStr, 'imflaw')
    imflaw = imageBlockProcess(im, [n_varblock n_varblock], ...
        o_varblock, @flawFun);
else
    imflaw = cacheStr.imflaw;
end

imflaw = imflaw - min(imflaw(:)); % Need non-negative cost map. Sometimes
                                  % round-off leaves us with slightly
                                  % negative numbers.
imflaw = imflaw / max(imflaw(:)); % Normalize for convenience


% Build the adjacency matrix
if or(~isfield(cacheStr, 'A'), ~isfield(cacheStr, 'ii') | ...
        ~isfield(cacheStr, 'bmin'))
    [A ii bmin] = mkAdjMat(imflaw, n_pathblock);
else
    A = cacheStr.A;
    ii = cacheStr.ii;
    bmin = cacheStr.bmin;
end

% Build a selection vector for only the border tiles
iborder = borderInOrder(ii);

% Calculate the cost map and predecessors for the border tiles
[d p] = calculateCuts(A, iborder, ii, b_pathblock);

% outPath = [];
% Reconstitute the minimum path

% figure(1); imagesc(imflaw); axis image; caxis([0 1]); hold on;

outPath = getPath(d, p, bmin, iborder);


% Cache outputs
if nargout > 1
    cacheStr.imflaw = imflaw;
    cacheStr.A = A;
    cacheStr.ii = ii;
    cacheStr.bmin = bmin;
    cacheStr.iborder = iborder;
    cacheStr.d = d;
    cacheStr.p = p;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n_varblock, o_varblock, n_pathblock, b_pathblock] = ...
    parseParms(cacheStr)
n_varblock = cacheStr.n_varblock;
o_varblock = cacheStr.o_varblock;
n_pathblock = cacheStr.n_pathblock;
b_pathblock = cacheStr.b_pathblock;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function op = getPath(d, p, bmin, iborder)
dmin = min(d(:));
[rm cm] = find(d == dmin, 1);
rmi = iborder(rm);
%cmi = iborder(cm);

p_coarse = [];
pr0 = rmi;
while pr0 > 0
    p_coarse = [p_coarse; pr0]; %#ok<AGROW>
    pr0 = p(pr0, cm);
end

[rmesh_bmin cmesh_bmin] = meshgrid(1:size(bmin,1), 1:size(bmin,2));
rc_coarse = cat(2, rmesh_bmin(p_coarse), cmesh_bmin(p_coarse));

op = cell(size(rc_coarse, 1) - 1, 1);

for k = 1:(size(rc_coarse, 1) - 1)
    op{k} = getFinePath(bmin, rc_coarse(k + 1,:), rc_coarse(k,:));
end

op = cat(1, op{:});

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [imrc] = getFinePath(bmin, rc1, rc2)
dv = rc2 - rc1;

blocksz = size(bmin(1).im);

% Reconstitute the sub-image, and the start and finish points within it
[im, rc_off1, rc_off2] = getcatim(bmin, rc1, dv);

[~, path] = getMinImagePath(im, rc_off1, rc_off2);

[rmesh cmesh] = meshgrid(1:size(im,2), 1:size(im,1));

imr = rmesh(path(:));
imc = cmesh(path(:));

rcoff = min(cat(1, rc1, rc2) - 1) .* blocksz;

imrc = cat(2, imr + rcoff(2), imc + rcoff(1));

end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d p] = calculateCuts(A, iborder, ii, b_pathblock)

d = zeros(size(A,1), numel(iborder));
p = d;
for k = 1:numel(iborder)
       [d(:,k), p(:,k)] = thirdpartyfunction('shortest_paths', A, iborder(k));
end

[rbo cbo] = meshgrid(1:numel(iborder), 1:numel(iborder));
blackout = abs(rbo - cbo) <= b_pathblock;

[rd cd] = meshgrid(1:size(ii,1), 1:size(ii,2));
rd = rd(iborder);
cd = cd(iborder);
dnorm = dist2(cat(2, rd, cd), cat(2, rd, cd));

d = d(iborder,:) ./ sqrt(dnorm);
d(blackout) = inf;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iborder = borderInOrder(ii)
iborder1 = ii(1:end,1);
iborder2 = ii(end,2:(end-1));
iborder3 = ii(end:-1:2,end);
iborder4 = ii(1,end:-1:2);

iborder = cat(1, iborder1(:), iborder2(:), iborder3(:), iborder4(:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A ii bmin] = mkAdjMat(imf, n)

bmin = imageBlockProcess(imf, [n n], 0, @blockmin, blockmin(0));

ii = 1:numel(bmin(:));
ii = reshape(ii, size(bmin));

% Find 8-connected graph weights
dv = [-1 -1; -1 0; -1 1; 0 1; 1 1; 1 0; 1 -1; 0 -1];

is = zeros([size(ii), size(dv, 1)]);
js = is;
as = is;

for r = 1:size(ii,1)
    for c = 1:size(ii,2)
        parfor d = 1:size(dv, 1)
            rd = r + dv(d,1);
            cd = c + dv(d,2); %#ok<PFBNS>
            
            if rcok(rd,cd,size(bmin))
                is(r, c, d) = ii(r, c); %#ok<PFBNS>
                js(r, c, d) = ii(rd, cd);
                [as(r, c, d) p] = getCost(bmin, [r c], [rd cd]);                
                as(r, c, d) = as(r, c, d) / numel(p);
            end
        end
        
        bmin(r, c).r = r;
        bmin(r, c).c = c;
    end
end

sel = is(:) > 0;

A = sparse(is(sel), js(sel), as(sel), numel(ii), numel(ii));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = rcok(r, c, sz)
rok = and(r > 0, r <= sz(1));
cok = and(c > 0, c <= sz(2));
ok = and(rok, cok);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cost path] = getCost(bmin, rc1, rc2)
dv = rc2 - rc1;

% Reconstitute the sub-image, and the start and finish points within it
[im, rc_off1, rc_off2] = getcatim(bmin, rc1, dv);

[cost path] = getMinImagePath(im, rc_off1, rc_off2);

end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [im, rc_o1, rc_o2] = getcatim(bmin, rc, dv)
o1 = [0 0];
o2 = [0 0];
rcd = rc + dv;
% two lines replace a whole bunch of if conditions. Maybe this is faster?
% if dv(1) < 0, br is [rc(1) - 1, rc(1)]
% if dv(1) > 0, br is [rc(1), rc(1) + 1]
% if dv(1) == 0, br is just rc(1).
% Similarly for dv(2) and bc.
br = (rc(1) - (dv(1) < 0)):(rc(1) + (dv(1) > 0));
bc = (rc(2) - (dv(2) < 0)):(rc(2) + (dv(2) > 0));

sub_bmin = bmin(br, bc);

% x-direction
if dv(2) < 0
    o1(2) = size(sub_bmin(1).im, 2);
else
    o2(2) = dv(2) * size(sub_bmin(1).im, 2);
end

% y
if dv(1) < 0
    o1(1) = size(sub_bmin(1).im, 1);
else
    o2(1) = dv(1) * size(sub_bmin(1).im, 1);
end
    
% rc in row-column, aka y-x
rc_o1 = [bmin(rc(1), rc(2)).y, bmin(rc(1), rc(2)).x] + o1;
rc_o2 = [bmin(rcd(1), rcd(2)).y, bmin(rcd(1), rcd(2)).x] + o2;



if size(sub_bmin,1) > 1
    im = [sub_bmin(1,:).im; sub_bmin(2,:).im];
else
    im = [sub_bmin.im];
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = blockmin(subim)
m = min(subim(:));

[o.y, o.x] = find(subim == m, 1);

o.im = subim;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = flawFun(subim)

d1 = diff(subim(:,2:end),1,1);
d2 = diff(subim(2:end,:),1,2);

d = d1.^2 + d2.^2;

dsamp = d(d <= ptile(d(:), .25));

m = sqrt(var(dsamp));

end
