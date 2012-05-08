function xyd = polynomialSDistortionApproximation(xy, k, o)

x = xy(:,1);
y = xy(:,2);

xd = x + y * k .* (3 / 8 + 3 / 4 * (x.^2 + y.^2) - ...
    1 / 8 * (x.^4 + y.^4) - 1 / 4 * x.^2 .* y.^2);
yd = y - x * k .* (3 / 8 + 3 / 4 * (x.^2 + y.^2) - ...
    1 / 8 * (x.^4 + y.^4) - 1 / 4 * x.^2 .* y.^2);

if o >= 2
    k2 = k.^2;
    xd = xd - x * k2 / 2 .* (x.^2 + y.^2);
    yd = yd - y * k2 / 2 .* (x.^2 + y.^2);
    
    if o >= 3
        k3 = k.^3;
        xd = xd + y * k3 .* (-3 / 48 * (x.^2 + y.^2) - 1 / 8 * (x.^4 + y.^4) -...
            1/4 * x.^2 .* y.^2);
        yd = yd - x * k3 .* (-3 / 48 * (x.^2 + y.^2) - 1 / 8 * (x.^4 + y.^4) -...
            1/4 * x.^2 .* y.^2);
    end
end

xyd = cat(2, xd, yd);
