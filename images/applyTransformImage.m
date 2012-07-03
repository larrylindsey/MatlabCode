function [im x y] = applyTransformImage(im, trans, x, y, interp)

if isempty(trans.inv)
    warning('applyTransformImage:Inverse', ...
        'Trans has no inverse, Populating it automatically');
    trans = fitInverseTransform(trans);
end

tr = maketform('custom', 2, 2, ...
    @(x, tr) applyTransform(x, tr.tdata), ...
    @(x, tr) applyTransform(x, tr.tdata.inv), trans);

if isfield(trans, 'Image') && ~isempty(trans.Image)
    doflip = true;    
    im = flipdim(im,1);
    disp('Got Reconstruct-style transform. Transforming image upside-down');
else 
    doflip = false;
end

args = {im, tr};

if nargin > 4
        args = [args {interp}];
end

if isfield(trans, 'data') && isfield(trans.data, 'u') && ...
        isfield(trans.data, 'v')
    args = [args, {'UData', trans.data.u, 'VData', trans.data.v}];   
end

if nargin > 2 && ~isempty(x) && ~isempty(y)
    
    args = [args, {'XData', x}];
    if nargin > 3
        args = [args, {'YData', y}];
    end
end

[im x y] = imtransform(args{:});

if doflip
    im = flipdim(im,1);
end

end
