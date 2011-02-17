function outgrp = groupByCycle(map)

outgrp = cell(size(map));

map = map(1:numel(map));
selmat = 1:numel(map);

donev = false(size(map));

nextmap = map;

while any(~donev)
    
    n = size(selmat, 1);
    
    mapBlock = repmat(nextmap, [n 1]);
    
    donev = donev | (sum(mapBlock == selmat, 1) > 0);
    
    selmat = cat(1, selmat, nextmap);

    nextmap = nextmap(map);
end

selmat = selmat(1:(end - 1), :);
grpsz = zeros(1, size(selmat, 2));

for i_t = 1:size(selmat, 2)
    outgrp{i_t} = unique(selmat(:,i_t));
    grpsz(i_t) = numel(outgrp{i_t});
end

[junk isort] = sort(grpsz, 'descend');

outgrp = {outgrp{isort}};

testv = 1:numel(grpsz);
ikeep = [];

while ~isempty(testv)
    i_t = testv(1);
    ikeep = [ikeep i_t];
    irm = [];

    for i_test = 1:numel(testv)
        i_i = testv(i_test);
        if any(intersect(outgrp{i_t}, outgrp{i_i}))
            outgrp{i_t} = unique([outgrp{i_t}; outgrp{i_i}]);
            irm = [irm i_test];
        end        
    end
    testv([1 irm]) = [];

end

outgrp = {outgrp{ikeep}};
end