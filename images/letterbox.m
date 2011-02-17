function imout = letterbox(im, sz)

sz = sz(1:2);
sz(3) = size(im, 3);

castfunc = str2func(class(im));

imout = castfunc(zeros(sz));

rows = size(im, 1);
cols = size(im, 2);

dr = sz(1) - rows;
dc = sz(2) - cols;

ri = (1:rows) + floor(dr / 2);
ci = (1:cols) + floor(dc / 2);

imout(ri, ci, :) = im;

end
