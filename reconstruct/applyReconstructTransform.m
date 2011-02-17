function [imout x y] = applyReconstructTransform(im, trinfo, x, y)

if size(im, 3) > 1
    im = rgb2gray(im);
end

im = im2double(im);

im = flipud(im);

if nargin > 2
    imout = imtransform(im, trinfo.tr, 'UData', trinfo.tr.tdata.data.u, ...
        'VData', trinfo.tr.tdata.data.v, 'XData', x, 'YData', y);
else
    [imout x y] = imtransform(im, trinfo.tr, 'UData', trinfo.tr.tdata.data.u, ...
        'VData', trinfo.tr.tdata.data.v);
end

imout = flipud(imout);

