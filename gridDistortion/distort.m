function [xout yout] = distort(x, y, K, x0,  y0)

d2 = (x - x0).^2 + (y - y0).^2;

xout = x + (x - x0) .* (K * d2);
yout = y + (y - y0) .* (K * d2);

end