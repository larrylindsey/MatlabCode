function m = modi(x, y)
if numel(y) > 1
    y = numel(y);
end

m = mod(x, y);

m(m == 0) = y;