function [uni extra] = makeFenceSample(h, t)
if nargin < 2
    t = 0;
end

L = bwlabel(h > t);

uni = zeros(max(L(:)), 2);
extra = zeros(numel(L), 3);
tt = zeros(numel(L), 1);
ll = zeros(size(tt));
hh = zeros(size(tt));

tt(:) = 1:numel(L);
ll(:) = L;
hh(:) = h;

extra(:, 1) = tt;
extra(:, 2) = ll;
extra(:, 3) = hh;

for i_l = 1:max(L(:))
    sel = logical(L == i_l);
    c = round(mean(tt(sel)));
    uni(i_l, 1) = c;
    uni(i_l, 2) = i_l;
end

end