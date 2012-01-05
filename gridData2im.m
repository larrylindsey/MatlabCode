function im = gridData2im(rc_grid, data)
% im = gridData2im(rc_grid, data)
%   rc_grid [n d] -  where n is the number of vectors, and d is the
%       dimesionality, which should probably be 2.
%   data [n] - the data to convert to image
%
%   im - an image such that im(rc_grid(i, 1), rc_grid(i, 2)) = data(i). For
%       points x, y where there is no i such that [x y] == rc_grid(i,:),
%       im(x, y) = nan.


[n d] = size(rc_grid);

rc_grid = rc_grid - repmat(min(rc_grid, [], 1), [n 1]) + 1;

sz = max(rc_grid, [], 1);

im = nan(sz);

im(sub2ind(sz, rc_grid(:,1), rc_grid(:,2))) = data;

end
