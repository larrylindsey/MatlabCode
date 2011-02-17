function Aout = taylorMatrix(x,y)
x1 = x(:);
y1 = y(:);

%x1, y1 in [n 1]

x2 = x1.*x1;
xy = x1.*y1;
y2 = y1.*y1;

x3 = x2.*x1;
x2y = x2.*y1;
xy2 = x1.*y2;
y3 = y2.*y1;

o = ones(size(x1));

Aout = cat(2, x3, x2y, xy2, y3, x2, xy, y2, x1, y1, o);