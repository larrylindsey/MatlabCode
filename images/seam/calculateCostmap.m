function costmap = calculateCostmap(energy)
%costmap = calculateCostmap(energy)
%
%Calculates the costmap of the given energy map.
%
% energy - the energy map for which a cost is to be calculated.
%
% costmap - a costmap for columnar-seam removal.  The value at position
%           (r,c) indicates the energy cost of removing the minimal seam
%           that connects the 1st row to (r,c).
sz = size(energy);

%Pre-allocate the costmap
costmap = zeros(sz);

%The first row is equal to the first row of the energy map.
costmap(1,:) = energy(1,:);

%Create row vectors where we'll store a row of costmap values, left shifted
%by one pixel and right shifted by one pixel.
%In the left shifted row, the rightmost value will be set as infinite, so
%that it won't affect the output of min.  Vice-versa for hte right shifted
%row.
costmap_row_ls = inf(1, sz(2));
costmap_row_rs = costmap_row_ls;

for r = 2:sz(1)    
    %Compute the shifted rows.
    costmap_row_rs(2:end) = costmap(r-1, 1:(end-1));
    costmap_row_ls(1:(end-1)) = costmap(r-1, 2:end);

    %Calculate the cost at row r.
    costmap(r, :) = energy(r,:) + min(...
        cat(1, costmap_row_ls, costmap(r-1,:), costmap_row_rs), ...
        [], 1);
end

