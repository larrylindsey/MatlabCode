function out = diffinterp1(in, xx, dim)

if nargin < 3
    sz = size(in);
    dim = find(sz > 1, 1, 'first');
    
    if isempty(dim)
        dim = 1;
    end
end

if nargin < 2 || isempty(xx)
    xx = zeros(size(in, dim), 1);
    xx(:) = 1:numel(xx);
end

xxndsz0 = ones(1, ndims(in));
xxndsz0(dim) = numel(xx);
xxnd = zeros(xxndsz0);
xxnd(:) = xx;

xxnd_repmatsz = size(in);
xxnd_repmatsz(dim) = 1;

xxnd = repmat(xxnd, xxnd_repmatsz);

inD = diff(in, 1, dim) ./ diff(xxnd, 1, dim);

dimperm = 1:ndims(inD);
dimperm(1) = dim;
dimperm(dim) = 1;

inD = permute(inD, dimperm);

dx = (xx(2:end) + xx(1:(end-1))) / 2;
out = interp1(dx, inD, xx, 'cubic', 'extrap');

out = ipermute(out, dimperm);