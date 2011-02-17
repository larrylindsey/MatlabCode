function [tri pts] = levelSetMesh(sets, z)

superset = cat(1, sets{:});

sel = cell(size(sets));
tri = [];
pts = zeros(size(superset, 1), 3);

n = 0;
for i_s = 1:numel(sets)
    m = size(sets{i_s}, 1);
    sel{i_s} = n + (1:m);
    
    pts(sel{i_s}, 1:2) = sets{i_s};
    pts(sel{i_s}, 3) = z(i_s);
    
    n = n + m;
end

for i_s = 2:numel(sets)
    map = mapLevelSetPair(sets{i_s - 1}, sets{i_s});
    remap = map;
    remap(:,1) = sel{i_s - 1}(map(:,1));
    remap(:,2) = sel{i_s}(map(:,2));
    preTri = mapToTriangles(remap);
    tri = cat(1, tri, preTri);
end

end

function tri = mapToTriangles(map)

n = size(map, 1);

tri = zeros(n, 3);

for i_m = 1:n
    ifwd = modi(i_m + 1, n);
    quad = map([i_m ifwd], :);
    newtri = unique(quad(:));
    if numel(newtri) > 3
        keyboard;
    end
    tri(i_m, :) = newtri;
end

end