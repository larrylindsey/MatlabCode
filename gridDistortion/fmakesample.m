function s = fmakesample(L)
v = 1:numel(L);
s = zeros(max(L), 2);
for i_l = 1:max(L)
    lsel = logical(L == i_l);
    s(i_l, 1) = mean(v(lsel));
    s(i_l, 2) = i_l;
end