function [dmap] = makeDistortionMap(trans, n, lim)

if nargin < 3
    lim = [0 1];
end
% 
% [R C] = meshgrid(linspace(lim(1), lim(2), 16), linspace(lim(1), lim(2), 16));
% RC = cat(2, R(:), C(:));
% RCtr = doTransform(RC, trans);
% 
% % ATrans.order = 1;
% % ATrans.T = calculateTransMat(RCtr, RC, ...
% %     ATrans.order, type);
% % ATrans.Tinv = calculateTransMat(RC, RCtr, ...
% %     ATrans.order, type);
% % ATrans.type = type;
% 
% 
% 
% %RCtr = doTransform(RCtr, ATrans);
% [RCtrLA sh r] = jankyLinearAlign(RCtr, RC, .1, .1);

[R C] = meshgrid(linspace(lim(1), lim(2), n), linspace(lim(1), lim(2), n));
RC = cat(2, R(:), C(:));
RCtr = doTransform(RC, trans);

RCtr = jankyLinearAlign(RCtr, RC);
% RCtr = RCtr + repmat(sh, [size(RCtr,1), 1]);
% zs = repmat(mean(RC,1), [size(RC,1) 1]);

% RCtr = RCtr - zs;
% RCtr = RCtr * r;
% RCtr = RCtr + zs;

dmap = reshape(sqrt(sum((RC - RCtr).^2, 2)), size(R));