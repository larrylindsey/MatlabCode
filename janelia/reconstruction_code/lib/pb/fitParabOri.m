function pbori2 = fitParabOri( pbori, thetalist, r )
% cut out of pbBGTG.m
%  Input
%       pbori -  y by x by orient
%       thetalist - 
%       r - savgol parameter
%   Output
%       pbori2 -    after parabolic sharpening,  y by x by orient

[h,w,norient] = size(pbori);
thetalist = pi*( 0:(norient-1) )/norient;
%[unused,maxo] = max(pbori,[],3);
pbori2 = zeros(h,w,norient);
%r = 2.5;
for i = 1:norient,
  %mask = (maxo == i);
  pbori2(:,:,i) = fitparab(pbori(:,:,i),r,r,thetalist(i));
end
pbori2 = max(0,min(1,pbori2));
