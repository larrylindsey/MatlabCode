function tform = xf_2_tform(xf, p, q)
% tform = xf_2_tform(xf, p, q)
% Generate imtransform compatible tform from a IMOD generated .xf file.
%
% Inputs:
%   xf            [1x6] vector of xf alignment transformation from IMOD
%   p             height of image.
%   q             width of image.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  05082008  Code borrowed from xf2tl.m [Dmitri Chklovskii].
%

tm(1:2,1:2)=[xf(1) -xf(3);-xf(2) xf(4)];
tm(3,1)=(1-xf(1))*q/2+xf(2)*p/2+xf(5);
tm(3,2)=(1-xf(4))*p/2+xf(3)*q/2-xf(6);
tform = maketform('affine',tm(:,:));

return
end
