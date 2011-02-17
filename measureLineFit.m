function e = measureLineFit(b, pts)
b = b(:);
c = pts * b - 1;
e = sqrt(mean(c.^2));