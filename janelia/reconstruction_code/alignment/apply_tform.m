function image_new = apply_tform(image, tg, margin, interpolation)
% image_new = apply_tform(image, tg, margin)
% Apply transform, tg, on an image keeping a margin around the border
%
% Inputs:
%   image           assumed to have single plane.
%   tg              TFORM compatible with imtransform.
%   margin          margin to be kept after transformed image.
%   interpolation   type of interpolation to be used during transformation.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  05082008  Code borrowed from newstack.m [Dmitri Chklovskii].
%

image_new = imtransform(image, tg, interpolation, 'XData', ...
  [margin.left margin.right],'YData',...
  [margin.top margin.bottom], 'FillValues', 0);

return
end
