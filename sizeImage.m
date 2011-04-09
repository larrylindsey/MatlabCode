function im = sizeImage(im, sz)


imsz = size(im);

if numel(imsz) == numel(sz) && all(imsz == sz)
    return;
end

if numel(imsz) > numel(sz)
    for d = (numel(sz) + 1):numel(imsz)
        im = mean(im, d);
    end
    imsz = size(im);
elseif numel(imsz) < numel(sz)
    repsz = sz;
    repsz(1:numel(imsz)) = 1;
    im = repmat(im, repsz);
    imsz = size(im);
end

tvar = im(1);
tvar(1) = 0;

for d = 1:numel(imsz)
    if imsz(d) < sz(d)
        padsize = sz(d) - imsz(d);
        headPadSz = imsz;
        tailPadSz = imsz;
        headPadSz(d) = floor(padsize / 2);
        tailPadSz(d) = ceil(padsize / 2);
        
        headPad = repmat(tvar, headPadSz);
        tailPad = repmat(tvar, tailPadSz);
        
        im = cat(d, headPad, im, tailPad);
    elseif imsz(d) > sz(d)
        cropsize = imsz(d) - sz(d);
        headIndex = ceil(cropsize / 2);
        tailIndex = headIndex - 1 + sz(d);
        selstr = makeSelectString(numel(imsz), d, headIndex, tailIndex);
        im = eval(sprintf('im(%s);', selstr));
    end
    imsz = size(im);
end

end

function selstr = makeSelectString(nd, d, llim, ulim)
selstr = '';
for id = 1:nd
    if id == d
        selstr = sprintf('%s,%d:%d', selstr, llim, ulim);
    else
        selstr = sprintf('%s,:', selstr);
    end
end

selstr = selstr(2:end);
end