function iout = linearIndex(sz, dim)
if size(dim, 2) ~= size(sz, 2)
    error('Size and Indices do not match');
end

dim = dim - 1;
nd = size(dim, 1);

szDot = cumprod([1 sz(1:(end - 1))]);

szDot = repmat(szDot, [nd, 1]);

iout = 1 + sum(szDot .* dim, 2);
