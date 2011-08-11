function seam = calculate_seam(energy, costmap)

sz = size(energy);

if nargin < 2
    costmap = calculateCostmap(energy);
end

totalcost = costmap(end,:);
seamc = zeros(1, sz(1));


[junk seamc(end)] = min(totalcost);

for r = (sz(1) - 1):-1:1
    seamc(r) = next_col(seamc(r+1), r, costmap);
end

seam = seamc;

end

function r_out = next_col(last_col, curr_row, costmap)
testcol_selector = (-1:1) + last_col;
testcol_selector = testcol_selector(logical(testcol_selector > 0));
testcol_selector = testcol_selector(logical(testcol_selector ...
    <= size(costmap, 2)));

[junk i_testcol] = min(costmap(curr_row, testcol_selector));
r_out = testcol_selector(i_testcol);
end