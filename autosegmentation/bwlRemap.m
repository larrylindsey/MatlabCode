function bwl3d = bwlRemap(bwl3d, eqmap)

if size(eqmap, 2) > 1
    vec = 1:max(bwl3d(:));
    for ii = 1:size(eqmap,1)
        vec(eqmap(ii,2)) = eqmap(ii,1);
    end
end

vec = [1 vec + 1];

bwl3d = vec(bwl3d + 1) - 1;


