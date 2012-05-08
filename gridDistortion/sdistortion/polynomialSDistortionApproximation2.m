function xyd = polynomialSDistortionApproximation2(xy, k, o)

x = xy(:,1);
y = xy(:,2);

xd = x + y * k .* (x.^2 + y.^2);
yd = y - x * k .* (x.^2 + y.^2);

if o >= 2
    k2 = k.^2;
    xd = xd - x * k2 .* ( (x.^4 + y.^4) / 2 + x.^2 .* y.^2 );
    yd = yd - y * k2 .* ( (x.^4 + y.^4) / 2 + x.^2 .* y.^2 );
    
end

xyd = cat(2, xd, yd);
