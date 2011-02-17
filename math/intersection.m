function [pt inBound] = intersection(line1, line2)

if isstruct(line1)
    x1 = line1.point1';
    x2 = line1.point2';
else
    x1 = line1(:,1);
    x2 = line1(:,2);
end

if isstruct(line2)
    y1 = line2.point1';
    y2 = line2.point2';
else
    y1 = line2(:,1);
    y2 = line2(:,2);
end

xhat = x2 - x1;
yhat = y2 - y1;

A = [xhat yhat];
b = x1 - y1;

if det(A) == 0
    pt = nan(2, 1);
    inBound = false;
    return;
end

v = inv(A) * b;

pt = v(2) * yhat + y1;

inBound = logical(v(1) < 0 && v(1) > -1 && v(2) > 0 && v(2) < 1);