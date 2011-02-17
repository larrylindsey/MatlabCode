function newslice = splitSlice(slice, newslice, n, f)

if nargin < 4
    f = 1;
    if nargin < 3
        n = 4;
    end
end

if nargin < 2 || isempty(newslice)
    newslice = uint8(zeros(size(slice)));
    istart = 1;
else
    istart = max(newslice(:)) + 1;
end

if isnumeric(n)
    s = strel('disk', n);
else
    s = n;
end

sliceErode = imerode(slice, s);
sliceEdge = mkEdge(slice);


nEdge = sum(sliceEdge(:));
nErode = sum(sliceErode(:));

factor = (nErode / nEdge / f);
factor = max(factor, 1);

erodeSel = find(sliceErode(:));
erodeDec = erodeSel;
erodeDec(round(1:factor:end)) = [];

sliceErode(erodeDec) = false;

newslice(sliceErode) = istart;
newslice(sliceEdge) = istart + 1;

end

function edge = mkEdge(slice)
r = size(slice, 1);
c = size(slice, 2);

Hposedge = diff(slice, 1, 1) == 1;
Hposedge = cat(1, false(1, c), Hposedge);

Hnegedge = diff(slice, 1, 1) == -1;
Hnegedge = cat(1, Hnegedge, false(1, c));

Vposedge = diff(slice,1 , 2) == 1;
Vposedge = cat(2, false(r, 1), Vposedge);

Vnegedge = diff(slice, 1, 2) == -1;
Vnegedge = cat(2, Vnegedge, false(r, 1));

edge = Hposedge | Hnegedge | Vposedge | Vnegedge;

end