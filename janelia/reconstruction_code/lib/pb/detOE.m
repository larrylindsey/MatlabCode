function [oe,thetalist] = detOE(im,sigma,elong,norient,deriv)
% modelled after:
% function [bg,tg,theta] = detBGTG(im,radius,norient)
%
% No Hilbert transf!
% No squaring!

if nargin<5, deriv=2; end

[h,w,unused] = size(im);
% filter bank (from fbCreate.m)
support = 3;
for orient = 1:norient,
    theta = (orient-1)/norient * pi;
    thetalist(orient) = theta;
    fb{orient} = oeFilter(sigma*[elong 1],support,theta, deriv); % 2 means deriv
    %no Hilbert fb{2*orient,scale} = oeFilter(sigma*[elong 1],support,theta,2,1);
end

% run
fim = fbRun( fb, im );
oe = zeros( h, w ); 
for o=1:numel(fim),
    oe(:,:,o) = fim{o};
end
