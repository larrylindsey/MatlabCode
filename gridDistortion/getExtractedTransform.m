function [sqtrans sq_rc_matcha trans rc_matcha] = ...
    getExtractedTransform(rc_found, rc_grid, grid_model, order, type, u, v)

if nargin < 7
    u = [min(rc_found(:,1)) max(rc_found(:,1))];
    v = [min(rc_found(:,2)) max(rc_found(:,2))];
    
    if nargin < 5
        type = @legendreMat;
    end
end

nansel = sum(isnan(rc_grid), 2) > 0;
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


sqtrans.order = order;
sqtrans.T = calculateTransMat(xy_found, sq_xy_matcha, order, type);
sqtrans.Tinv = calculateTransMat(sq_xy_matcha, xy_found, order, type);
sqtrans.iDim = 2;
sqtrans.oDim = 2;
sqtrans.type = type;
sqtrans.data.u = u;
sqtrans.data.v = v;
sqtrans.data.n = size(rc_found, 1);


trans.order = order;
trans.T = calculateTransMat(xy_found, xy_matcha, order);
trans.Tinv = calculateTransMat(xy_matcha, xy_found, order);
trans.iDim = 2;
trans.oDim = 2;
trans.type = type;
trans.data  = sqtrans.data;

end

function [mag angle] = polarModel(grid_model)
mag = sqrt(abs(det(grid_model)));

ngm = grid_model / mag;
v1 = ngm * [1 ; 0];
v2 = [0 1 ; -1 0] * ngm * [0 ; 1];

v = (v1 + v2) / 2;

angle = atan2(v(2), v(1));
end