function imbw = removeIslands(imbw, t, n)

ll = bwlabel(imbw, n);

for il = 1:max(ll(:))
    sel = ll == il;
    ct = sum(sel(:));
    if ct <= t
        imbw(sel) = false;
    end
end

end
