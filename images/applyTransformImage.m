function [im x y] = applyTransformImage(im, trans, x, y, interp)

tr = maketform('custom', trans.iDim, trans.oDim, @doCubicTransform, ...
    @doInverseCubicTransform, trans);

args = {im, tr};

if nargin > 4
        args = {args{:} interp};
end

if isfield(trans, 'data') && isfield(trans.data, 'u') && ...
        isfield(trans.data, 'v')
    args = {args{:}, 'UData', trans.data.u, 'VData', trans.data.v};   
end

if nargin > 2
    
    args = {args{:}, 'XData', x};
    if nargin > 3
        args = {args{:}, 'YData', y};
    end
end

[im x y] = imtransform(args{:});
