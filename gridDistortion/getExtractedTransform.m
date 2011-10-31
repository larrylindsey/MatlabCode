function output = getExtractedTransform(rc, cgrid, gm, control)
% function tr = getExtractedTransform(rc, cgrid, gm, control)

% return an example control struct
if nargin < 1
    tr.order = 3;
    tr.u = [-1 1];
    tr.v = [-1 1];
    tr.norm = [];
    tr.type = @legendreMat;
    tr.useRansac = true;
    output = tr;
    return;
end

if det(gm) < 0
    gm = gm([2 1], :);
    cgrid = cgrid(:,[2 1]);
end

% Do normalization, if specified.
if isfield(control, 'norm') && ~isempty(control.norm)
    norm = control.norm;
    if numel(norm) == 1
        norm = repmat(norm, [1 2]);
    end
    
    gm = 2 * gm ./ repmat(norm, [2 1]);
    
    norm = repmat(norm, [size(rc, 1), 1]) / 2;
    rc = (rc ./ norm) - 1;
    
end

% square the grid model
squareGM = getSquareGM(gm);

% Calculate the approximated true grid
squareRC = getSquareRC(squareGM, cgrid, rc);

data.u = control.u;
data.v = control.v;
data.n = round(sqrt(size(rc, 1)));

if isfield(control, 'useRansac') && control.useRansac
    ctrl_ransac = ransacRegressionTransform();
    ctrl_ransac.n = 24;
    ctrl_ransac.minInliers = 30;
    ctrl_ransac.maxIter = 100;
    [trsim sel] = ransacRegressionTransform(ctrl_ransac, ...
        squareRC, rc, control.order, control.type, data);%#ok
else
    sel = true(size(rc,1), 1);
end

% Reject outliers from the rc's
[traff affRC] = getTRAff(rc(sel,:), squareRC(sel,:), data, control);
squareRC = trsAlign(squareRC(sel,:), rc(sel,:));

trsim = regressionTransform(rc(sel,:), squareRC, control.order, ...
    control.type, data);

% traff - transform from rc to affine-aligned "square" grid
% trsim - transform from rc to similiarity-aligned square grid


rc_aff = doTransform(rc(sel,:), traff);
rc_sim = doTransform(rc(sel,:), trsim);

output.rc_square = squareRC;
output.rc_affine = affRC;
output.similarity.tr = trsim;
output.similarity.err =  rms(sum((rc(sel,:) - rc_sim).^2, 2));
output.similarity.rc = rc_sim;
output.affine.tr = traff;
output.affine.err =  rms(sum((rc(sel,:) - rc_aff).^2, 2));
output.affine.rc = rc_aff;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [traff aff_rc] = getTRAff(rc, squareRC, data, control)

tral = regressionTransform(squareRC, rc, 1, @legendreMat, data);
aff_rc = doTransform(squareRC, tral);

traff = regressionTransform(aff_rc, rc, control.order, control.type,...
    data);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sqgm = getSquareGM(gm)
v1 = gm(1,:);
v2 = gm(2,:);
v1n = norm(v1);
v2n = norm(v2);

cosups = (v1(1) * v2n + v2(2) * v1n) / (2 * v2n * v1n);
sinups = sqrt(1 - cosups*cosups);

sqgm = [cosups sinups; -sinups cosups];

if v1(2) < 0
    sqgm = sqgm';
end

gmnorm = (norm(gm(1,:)) + norm(gm(2,:))) / 2;

sqgm = sqgm * gmnorm;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function squareRC = getSquareRC(squareGM, cgrid, rc)

squareRC = cgrid * (squareGM');
squareRC = squareRC + ...
    repmat(mean(rc, 1) - mean(squareRC, 1), [size(squareRC, 1), 1]);

squareRC = trsAlign(squareRC, rc);
end