function output = extractedTransform(rc_found, rc_grid, grid_model, order)
% output = extractedTransform(rc_found, rc_grid, grid_model)
%
%

% Note: this function is a rewrite of getExtractedTransform


gmSQ = getSquareGM(grid_model);

% squareMatchGrid is a set of locations on a square grid,
% similarity-aligned to rc_found
squareMatchGrid = getSquareMatchGrid(rc_found, gmSQ, rc_grid);

% affineMatchGrid as above, but affine-aligned
affineMatchGrid = getAffineMatchGrid(rc_found, gmSQ, rc_grid);

square_str = repmat(struct, [1 numel(order)]);
affine_str = repmat(struct, [1 numel(order)]);


for ii = 1:numel(order)
    
    trSQ = regressionTransform(rc_found, squareMatchGrid, order(ii));
    trAff = regressionTransform(rc_found, affineMatchGrid, order(ii));
    
    rcSQMatch = doTransform(rc_found, trSQ);
    rcAffMatch = doTransform(rc_found, trAff);
        
    square_str(ii).rc_grid = squareMatchGrid;
    square_str(ii).tr = trSQ;
    square_str(ii).rc_match = rcSQMatch;
    square_str(ii).order = order(ii);
    
    affine_str(ii).rc_grid = affineMatchGrid;
    affine_str(ii).tr = trAff;
    affine_str(ii).rc_match = rcAffMatch;
    affine_str(ii).order = order(ii);
        
end

output.gm = grid_model;
output.gmSQ = gmSQ;
    
% output.square.rc_grid = squareMatchGrid;
% output.square.tr = trSQ;
% output.square.rc_match = rcSQMatch;
% 
% output.affine.rc_grid = affineMatchGrid;
% output.affine.tr = trAff;
% output.affine.rc_match = rcAffMatch;

output.square = square_str;
output.affine = affine_str;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gmSQ = getSquareGM(grid_model)

gmDet = det(grid_model);
d = sqrt(abs(gmDet));

gmSQ = eye(2) * d;

if gmDet < 0
    warning('extractedTransform:swappedAxes',...
        'Detected swapped axes in grid model');
    gmSQ = gmSQ([2 1],:);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function matchGrid = getSquareMatchGrid(rc_pt, gm, rc_grid)

grid0 = rc_grid * gm;
matchGrid = trsAlign(grid0, rc_pt);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function matchGrid = getAffineMatchGrid(rc_pt, gm, rc_grid)

grid0 = rc_grid * gm;
tr = regressionTransform(grid0, rc_pt, 1);
matchGrid = doTransform(grid0, tr);

end