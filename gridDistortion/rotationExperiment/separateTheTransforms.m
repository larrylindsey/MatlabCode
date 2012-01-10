function [trG trM trR tr00str tr90str] = separateTheTransforms(rot00, rot90,...
    gridRot)
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
global data GD_TOL dT K N_KTEST;
%constants
ORDER = 5;
GD_TOL = 1e-6;
dT = 1e-6;
N_KTEST = 24;
% K = 1e-1;
K = 16;
N = 32;

control = getExtractedTransform;
control.order = ORDER;
control.useRansac = false;

data.u = [-1 1];
data.v = [-1 1];
data.n = 32;
data.gamma = [];
data.weight = [];
data.order = ORDER;
data.type = control.type;


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
trMRM = regressionTransform(rc, rcMRM, control.order, control.type, data);
trGRG = regressionTransform(rc, rcGRG, control.order, control.type, data);

trG = invertTransStruct(estimateTransform(rc, trGRG, trR, control));
trM = estimateTransform(rc, trMRM, trR, control);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trT = estimateTransform(rc, trTRT, trR, control)
global data GD_TOL K;
kmag = K;
tr = regressionTransform(rc, rc, round(control.order / 2), control.type, data);

cerr = transformError(tr, trR, rc, trTRT, control, data);
dd = inf;

while dd > GD_TOL || dd < 0;
    lerr = cerr;
    grad = estimateGrad(tr, trR, rc, trTRT, control, data);
    ltr = tr;
    
    [k kmag] = estimateK(kmag, tr, grad, trR, rc, trTRT, control);
    tr = updateTR(tr, grad, k);
    
    cerr = transformError(tr, trR, rc, trTRT, control, data);
    fprintf('Last error\t%g\nCurr error\t%g\nTol\t\t%g\nk\t\t%g\n\n', ...
        lerr, cerr, GD_TOL, k);
    
    dd = lerr - cerr;
    
    if dd < 0 %&& k * 128 > 1
         kmag = kmag / 2;
        cerr = lerr;
        tr = ltr;
    end
end

fprintf('DONE!\n');

trT = tr;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [k kmag] = estimateK(kmag, tr, grad, trR, rc, trTRT, control)
global N_KTEST data;

pdata = data;
ktest = [0 logspace(-6, log(kmag)/log(10), N_KTEST-1)];
etest = zeros(size(ktest));

parfor ii = 1:numel(ktest)
    etest(ii) = transformError(updateTR(tr, grad, ktest(ii)), trR, rc, trTRT,...
        control, pdata);
end

[~, imin] = min(etest);

interpsel = imin + [-1 0 1];
if interpsel(1) < 1
    interpsel = 1:3;
end
if interpsel(3) > numel(etest)
    interpsel = numel(etest) + [-2 -1 0];
end

k = interpQuadraticExtremum(ktest(interpsel), etest(interpsel));

figure(1);
subplot(3,1,1);
plot(ktest, etest);
subplot(3,1,2);
plot(grad(:));
subplot(3, 1, 3);
plot(tr.T(:));
drawnow;
% if k <= ktest(round(numel(ktest)/4))
%     kmag = kmag / 2;
% end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grad = estimateGrad(tr, trR, rc, trTRT, control, data)
global dT;
cdt = dT;
grad = zeros(size(tr.T));
trTest = repmat(tr, [1 numel(grad)]);

parfor n = 1:numel(grad)
    trTest(n).T(n) = trTest(n).T(n) + cdt; %#ok<PFOUS>
    trTest(n) = populateTransInverse(trTest(n));
    perr = transformError(trTest(n), trR, rc, trTRT, control, data);
    
    trTest(n).T(n) = trTest(n).T(n) - cdt;
    trTest(n) = populateTransInverse(trTest(n));
    nerr = transformError(trTest(n), trR, rc, trTRT, control, data);
    
    grad(n) = (perr - nerr) / 2 / cdt;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tr = updateTR(tr, grad, k)
tr.T = tr.T - k * grad;
tr = populateTransInverse(tr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = transformError(tr, trR, rc, trTRT, control, data)

rct = doTransform(rc, tr);
rct = doTransform(rct, trR);
rct = doTransform(rct, invertTransStruct(tr));
trRCT = regressionTransform(rc, rct, control.order, control.type, data);

e = mean((sum((trTRT.T(:) - trRCT.T(:)).^2,2)));
%e = rms(trRCT.T(:) - trTRT.T(:));

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
