function [im_out sz_out] = hThenWResize(im, reduction, varargin)
[im_out sz] = reduceHeight(im, reduction(1), varargin{:});
[im_out sz2] = reduceWidth(im_out, reduction(2), varargin{:});
sz_out = cat(1, sz, sz2);
