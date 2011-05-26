function [trsim traff rc squareRC] = getExtractedTransform(rc, cgrid, gm, control)
% function tr = getExtractedTransform(rc, cgrid, gm, control)

% return an example control struct
if nargin < 1
    trsim.order = 3;
    trsim.u = [-1 1];
    trsim.v = [-1 1];
    trsim.norm = [];
    trsim.type = @legendreMat;
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
squareRC = cgrid * (squareGM');
squareRC = squareRC + ...
    repmat(mean(rc, 1) - mean(squareRC, 1), [size(squareRC, 1), 1]);

squareRC = trsAlign(squareRC, rc);

data.u = control.u;
data.v = control.v;
data.n = round(sqrt(size(rc, 1)));
%tr = regressionTransform(squareRC, rc, control.order, control.type, data);
[trsim sel] = ransacRegressionTransform(ransacRegressionTransform, ...
    squareRC, rc, control.order, control.type, data);%#ok

traff = getTRAff(rc, squareRC, sel, data, control);

squareRC = trsAlign(squareRC(sel,:), rc(sel,:));
trsim = regressionTransform(squareRC, rc(sel,:), control.order, control.type, ...
    data);

% % Fix the mag
% tr.T(2:3,:) = tr.T(2:3,:) / (sqrt(det(tr.T(2:3,:)) * 3 / 4 ));
% tr = populateTransInverse(tr);

end

function traff = getTRAff(rc, squareRC, sel, data, control)
rc = rc(sel,:);
squareRC = squareRC(sel,:);

tral = regressionTransform(squareRC, rc, 1, @legendreMat, data);
squareRC = doTransform(squareRC, tral);

traff = regressionTransform(squareRC, rc, control.order, control.type,...
    data);
end

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