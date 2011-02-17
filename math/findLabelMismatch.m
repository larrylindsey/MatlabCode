function [map1 map2] = findLabelMismatch(label1, label2)

map1 = {};
map2 = {};

h = waitbar(0, 'doing stuff');

m = max(label1(:));

for i_l = 1:m
    sel1 = logical(label1 == i_l);
    v2 = label2(sel1);
    v2 = v2(1);
    sel2 = logical(label2 == v2);
    if any(xor(sel1, sel2))
        map1 = cat(1, map1, {i_l});
        i21 = label2(sel2);
        i22 = label2(sel1);
        map2 = cat(1, map2, unique(cat(1, i21, i22)));
    end
    waitbar(i_l / m, h);
end

disp('Done.');

close(h);