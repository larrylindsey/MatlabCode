function [f a] = countContourXSec(contour)

sec = [contour.section];
maxSec = max(sec);

f = zeros(1, maxSec);
acell = cell(1, maxSec);

for s = 1:maxSec
    sel = find(sec == s);
    f(s) = numel(sel);
    acell{s} = zeros(1, f(s));
    for i_f = 1:f(s)
        c = contour(sel(i_f));
        acell{s}(i_f) = polyarea(c.transPoints(:,1), c.transPoints(:,2));
    end
end

a = zeros(max(f), maxSec);

for s = 1:maxSec
    ac = acell{s};
    for ia = 1:numel(ac)
        a(ia, s) = ac(ia);
    end
end