function [trG trM tr00str tr90str] = separateTheTransforms(rot00, rot90, gridRot)
% function tr = separateTheTransforms(rot00, rot90, gridRot)
% 
% Attempts to use two calibration images, one rotated and/or translated
% from the other, to separate the "calibration transform" from the
% "microscope transform"
%
% rot00 and rot90 are structs with fields rc_found, rc_grid, grid_model,
% and scale, as returned by extractDistortionGrid
%
% 
% rot00 is assumed to be representable by Tm o Tg, where Tm and Tg are
% functions from R2 to R2, Tm representing microscope distortion, and Tg
% representing sample-intrinsic distortion (ie, a map from a nominally
% ideal, perfectly square cross grating to the sample as it exists in
% reality).
%
% rot90 is assumed to be representable by Tm o Tr o Tg, where Tm and Tg are
% as above, and Tr is a rotation + translation.
global data GD_TOL dT K;
%constants
GD_TOL = 1e-6;
dT = 1e-6;
% K = 1e-1;
K = 1;
ORDER = 5;
N = 32;
data.u = [-1 1];
data.v = [-1 1];
data.n = 32;
data.gamma = [];
data.weight = [];

control = getExtractedTransform;
control.order = ORDER;
control.useRansac = false;

if nargin < 3
    gridRot = 1;
end

[tr00str tr90str rcFound00, rcFound90] =...
    calculateCompoundTransforms(rot00, rot90, gridRot, control);

trR = calculateR(rcFound00, rcFound90);

rc = gridRC(linspace(-1,1,N), linspace(-1,1,N));
rcMRM = sepX(rc, invertTransStruct(tr00str.similarity.tr),...
    tr90str.similarity.tr);
rcGRG = sepX(rc, tr90str.similarity.tr,...
    invertTransStruct(tr00str.similarity.tr));

trG = invertTransStruct(estimateTransform(rc, rcGRG, trR, control));
trM = estimateTransform(rc, rcMRM, trR, control);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trT = estimateTransform(rc, rcTRT, trR, control)
global data GD_TOL K;
k = K;
tr = regressionTransform(rc, rc, round(control.order / 2), control.type, data);

cerr = transformError(tr, trR, rc, rcTRT);
dd = inf;

while dd > GD_TOL || dd < 0;
    lerr = cerr;
    grad = estimateGrad(tr, trR, rc, rcTRT);
    ltr = tr;
    tr = updateTR(tr, grad, k);
    
    cerr = transformError(tr, trR, rc, rcTRT);
    fprintf('Last error\t%g\nCurr error\t%g\nTol\t\t%g\nk\t\t%g\n\n', ...
        lerr, cerr, GD_TOL, k);
    
    dd = lerr - cerr;
    
    if dd < 0 %&& k * 128 > 1
        k = k / 2;
        cerr = lerr;
        tr = ltr;
    end
end

fprintf('DONE!\n');

trT = tr;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grad = estimateGrad(tr, trR, rc, rcTRT)
global dT;
cdt = dT;
grad = zeros(size(tr.T));
trTest = repmat(tr, [1 numel(grad)]);

parfor n = 1:numel(grad)
    trTest(n).T(n) = trTest(n).T(n) + cdt; %#ok<PFOUS>
    trTest(n) = populateTransInverse(trTest(n));
    perr = transformError(trTest(n), trR, rc, rcTRT);
    
    trTest(n).T(n) = trTest(n).T(n) - cdt;
    trTest(n) = populateTransInverse(trTest(n));
    nerr = transformError(trTest(n), trR, rc, rcTRT);
    
    grad(n) = (perr - nerr) / 2 / cdt;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tr = updateTR(tr, grad, k)
tr.T = tr.T - k * grad;
tr = populateTransInverse(tr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = transformError(tr, trR, rc, rcTRT)
rct = doTransform(rc, tr);
rct = doTransform(rct, trR);
rct = doTransform(rct, invertTransStruct(tr));

e = mean((sum((rcTRT - rct).^2,2)));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rcTRT = sepX(rc, trA, trB)
rcA = doTransform(rc, trA);
rcTRT = doTransform(rcA, trB);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trR = calculateR(rcFound00, rcFound90)
%half assed alignment
global data;

rc00trs = trsAlign(rcFound00, rcFound90);
trR = regressionTransform(rcFound00, rc00trs, 1, @legendreMat, data);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tr00str tr90str rcFound00 rcFound90] = ...
    calculateCompoundTransforms(rot00, rot90, gridRot, control)

switch gridRot
    case -1
        GR = [0 1; -1 0];
    case 0
        GR = eye(2);
    case 1
        GR = [0 -1; 1 0];
    otherwise
        error('gridRot must be -1, 0, or 1');
end

% Calculate tr00 and tr90 from the discovered grid distortions
[sel00 sel90] = rcMatch(rot00.rc_grid, rot90.rc_grid * GR);

rcGrid00 = rot00.rc_grid(sel00,:);
rcGrid90 = rot90.rc_grid(sel90,:);

rcFound00 = rot00.rc_found(sel00,:);
rcFound90 = rot90.rc_found(sel90,:);

gm = rot00.grid_model;

tr00str = getExtractedTransform(rcFound00, rcGrid00, gm, control);
tr90str = getExtractedTransform(rcFound90, rcGrid90, gm, control);
end
