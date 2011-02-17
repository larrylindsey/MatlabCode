function [match dd ng] = gridEstimate(gridStruct, im)

fillgap_len = max(size(im)) / 16;
grid_size_threshold = 32;

[deg bw h t r peaks] = parseInput(gridStruct);


lines = houghlines(bw, t, r, peaks, 'FillGap', fillgap_len);

tlines = [lines.theta];
tlines = mod(tlines - deg + 45, 180) - 45;

rlines = [lines.rho];

ihoriz = logical(abs(tlines) < 3);
ivert = logical(abs(tlines - 90) < 3);

lines_horiz = rejectIntersectors(lines(ihoriz), size(im));
lines_vert = rejectIntersectors(lines(ivert), size(im));

ddh = diff(sort(rlines(ihoriz)));
ddv = diff(sort(rlines(ivert)));

ddh = ddh(logical(ddh ~= 0));
ddv = ddv(logical(ddv ~= 0));

dh = median(ddh(ddh > grid_size_threshold));
dv = median(ddv(ddv > grid_size_threshold));

dd = mean([dh dv]);
dd(2) = dd;

hnum = size(im, 2) / dh;
vnum = size(im, 1) / dv;

ng = [hnum vnum];


pts = zeros(2, numel(lines_vert), numel(lines_horiz));

parfor i_h = 1:numel(lines_horiz)
    ptsh = zeros(2, numel(lines_vert));
    for i_v = 1:numel(lines_vert)
        ptsh(:, i_v) = intersection(lines_vert(i_v), ...
            lines_horiz(i_h));
    end
    pts(:,:,i_h) = ptsh;
end

pts = reshape(pts, [2 numel(lines_horiz) * numel(lines_vert)]);
pts = removeRedundantPoints(pts, grid_size_threshold);

pts = round(pts);


gsz = ceil(dd);

igood = (pts(1,:) > gsz(2)) .* (pts(1,:) < size(im, 2) - gsz(2)) .* ...
    (pts(2,:) > gsz(1)) .* (pts(2,:) < size(im, 1) - gsz(1));
igood = logical(igood);

pts = pts(:, igood);

rsel = -gsz(1):gsz(1);
csel = -gsz(2):gsz(2);

castfunc = str2func(class(im));

patches = castfunc(zeros(gsz(1) * 2 + 1, gsz(2) * 2 + 1, size(pts, 2)));

for i_p = 1:size(pts, 2)
    patches(:, :, i_p) = im(pts(2, i_p) + rsel, ...
        pts(1, i_p) + csel);    
end

match = median(patches, 3);

%match = imrotate(match, deg, 'bilinear', 'crop');

%msz = size(match);
%m = floor(msz * .1 / 2);
%match = match(m(1):(msz(1) - m(1)), m(2):(msz(2) - m(2)));
end

function [deg bw h t r peaks sz] = parseInput(gridStruct)

deg = gridStruct.deg;
bw = gridStruct.bw;
h = gridStruct.h;
t = gridStruct.t;
r = gridStruct.r;
peaks = gridStruct.peaks;
sz = gridStruct.sz;
end

function l = lineLength(line)
l = sqrt(sum((line.point2 - line.point1).^2));
end

function lines = rejectIntersectors(lines, rcsz)
keepsel = true(1, numel(lines));

n = numel(lines);

px = [0 0 rcsz(2) rcsz(2)];
py = [0 rcsz(1) rcsz(1) 0];

for ii = 1:n

    jj = ii + 1;
    
    while keepsel(ii) && jj <=n
        pt = intersection(lines(ii), lines(jj));
        if inpolygon(pt(1), pt(2), px, py)
            if lineLength(lines(ii)) > lineLength(lines(jj))
                keepsel(jj) = false;
            else
                keepsel(ii) = false;
            end
        end
        jj = jj + 1;
    end
    
end

end

function pts = removeRedundantPoints(pts, t)
pts = pts';
n = size(pts, 1);
sel = 1:n;

tsqr = t * t;

isel = 1;
while isel < numel(sel)
    ii = sel(isel);
    
    loopsel = sel(sel > ii);    
    dmap = dist2(pts(ii,:), pts(loopsel, :));
    
    equivPtsSel = loopsel(dmap < tsqr);
    
    newpt = mean(pts([ii; equivPtsSel(:)], :), 1);
    pts(ii,:) = newpt;
    
    sel = setdiff(sel, equivPtsSel);
    
    isel = isel + 1;
end

pts = pts(sel,:);

pts = pts';
end