function [sq_rc_matcha rc_matcha sel] = ...
    getMatchGrid(rc_found, rc_grid, grid_model, type)

    
if nargin < 4
    type = @legendreMat;
end


nansel = sum(isnan(rc_grid), 2) > 0;
sel = not(nansel);
rc_grid = rc_grid(not(nansel),:);
xy_found = rc_found(not(nansel),[2 1]);
xy_grid = rc_grid(:,[2 1]);

[mag angle] = polarModel(grid_model);

sq_grid_model = (eye(2) * mag) * ...
    [cos(angle) sin(angle); -sin(angle) cos(angle)];

sq_xy_match = xy_grid * sq_grid_model;
xy_match = xy_grid * grid_model;

% Affine alignment
sqtrans_align.order = 1;
sqtrans_align.T = calculateTransMat(sq_xy_match, xy_found, ...
    sqtrans_align.order, type);
sqtrans_align.Tinv = calculateTransMat(xy_found, sq_xy_match, ...
    sqtrans_align.order, type);
sqtrans_align.type = type;

trans_align.order = sqtrans_align.order;
trans_align.T = calculateTransMat(xy_match, xy_found, ...
    sqtrans_align.order, type);
trans_align.Tinv = calculateTransMat(xy_found, xy_match, ...
    sqtrans_align.order, type);
trans_align.type = type;


sq_xy_matcha = doTransform(sq_xy_match, sqtrans_align);
xy_matcha = doTransform(xy_match, trans_align);

rc_matcha = xy_matcha(:,[2 1]);
sq_rc_matcha = sq_xy_matcha(:,[2 1]);

end


function [mag angle] = polarModel(grid_model)
mag = sqrt(abs(det(grid_model)));

ngm = grid_model / mag;
v1 = ngm * [1 ; 0];
v2 = [0 1 ; -1 0] * ngm * [0 ; 1];

v = (v1 + v2) / 2;

angle = atan2(v(2), v(1));
end