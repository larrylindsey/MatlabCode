function [im_vect mag] = vectorMapToImage(u, v, n)

szu = size(u);
szv = size(v);

if any(szu ~= szv)
    error('x, y, u, and v must all have the same size');
end

mag = sqrt(u.^2 + v.^2);
atans = atan2(v, u);

if nargin < 3 || isempty(n)
    n = max(mag(:));
    disp(n);
end

if numel(n) > 1
    n = max(mag(n(:)));
end

hues = wheelHue(atans);
sat = mag / n;
sat(sat > 1) = 1;
val = ones(size(mag));

hsv = cat(2, hues(:), sat(:), val(:));
rgb = hsv2rgb(hsv);

im_vect = reshape(rgb, [size(u) 3]);
