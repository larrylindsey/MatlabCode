function x0 = interpQuadraticExtremum(x, y)

x12_2 = x(1)^2 - x(2)^2;
x13_2 = x(1)^2 - x(3)^2;
x23_2 = x(2)^2 - x(3)^2;

x12 = x(1) - x(2);
x13 = x(1) - x(3);
x23 = x(2) - x(3);

xext_num = y(1) * x23_2 - y(2) * x13_2 + y(3) * x12_2;
xext_den = y(1) * x23 - y(2) * x13 + y(3) * x12;
x0 = .5 * xext_num / xext_den;

end
