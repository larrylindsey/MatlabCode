function [im_out sz_out] = wThenHResize(im, reduction, varargin)
[im_out sz] = reduceWidth(im, reduction(2), varargin{:});
[im_out sz2] = reduceHeight(im_out, reduction(1), varargin{:});
sz_out = cat(1, sz, sz2);
