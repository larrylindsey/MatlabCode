function pb2 = gentleNonmax( pbori, thetalist )
% cut out of makePB2.m
% nonmax suppression and max over orientations
%  Input
%       pbori -  oriented pb, y by x by orient
%       thetalist - obsolete
%   Output
%       pb2   -  oriented pb, y by x by orient, with nonmax values set to 0

[h,w,norient] = size(pbori);
thetalist = pi*( 0:(norient-1) )/norient;

pb2 = zeros(h,w,norient);
for i = 1:norient,
  pb2(:,:,i) = nonmax( pbori(:,:,i), thetalist(i) );
end
