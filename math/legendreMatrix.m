function Aout = legendreMatrix(x,y)
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

p0 = o / 2;

p10 = sqrt(3/4) * x1;
p01 = sqrt(3/4) * y1;

p20 = .75 * sqrt(5) * (x2 - 1/3);
p11 = 1.5 * xy;
p02 = .75 * sqrt(5) * (y2 - 1/3);

p30 = (5/4) * sqrt(7) * (x3 - (3/5) * x1);
p21 = (3/4) * sqrt(15) * (x2y - y1/3);
p12 = (3/4) * sqrt(15) * (xy2 - x1/3);
p03 = (5/4) * sqrt(7) * (y3 - (3/5) * y1);

Aout = cat(2, p30, p21, p12, p03, p20, p11, p02, p10, p01, p0);