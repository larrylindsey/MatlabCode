function alignedCompositeImages(cache, imfiles, outdir, preTr)

unix(['mkdir -p ' outdir]);

if nargin < 4
    preTr = [];
end


ncpu = max(matlabpool('size'), 1);
splitImfiles0 = cell(1, ncpu);
splitImfiles1 = cell(1, ncpu);
splitTr0 = cell(1, ncpu);
splitTr1 = cell(1, ncpu);
splitIndices = cell(1, ncpu);

trAlign = seriesAlignmentTransforms(cache, round(numel(cache.feat1) / 2), ...
    preTr);

[xbound ybound] = getBounds(trAlign);

for ii = 1:ncpu
    splitIndices{ii} = ii:ncpu:numel(cache.feat1);
    splitTr0{ii} = trAlign(splitIndices{ii});
    splitTr1{ii} = trAlign(splitIndices{ii} + 1);
    splitImfiles0{ii} = imfiles(splitIndices{ii});
    splitImfiles1{ii} = imfiles(splitIndices{ii} + 1);
end



parfor cpu = 1:ncpu
    tr1 = splitTr0{cpu};
    tr2 = splitTr1{cpu};
    imfiles0 = splitImfiles0{cpu};
    imfiles1 = splitImfiles1{cpu};
    
    for ii = 1:numel(tr1)
        im1 = getIm(imfiles0{ii});
        im2 = getIm(imfiles1{ii});
        
        im1 = applyTransformImage(im1, tr1(ii), xbound, ybound);
        im2 = applyTransformImage(im2, tr2(ii), xbound, ybound);
        
        im1 = clipIm(im1);
        im2 = clipIm(im2);
        
        imwrite(compositeImage(im1, im2)), ...
            sprintf('%s/th %03d aligned %03d.png', outdir, ...
            indices(ii), indices(ii) + 1);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xbound ybound] = getBounds(trAlign)

xbound = [inf, -inf];
ybound = xbound;
n = numel(trAlign);

for ii = 1:n
    u = trAlign(ii).data.u;
    v = trAlign(ii).data.v;
    UV = gridRC(u,v);
    XY = applyTransform(UV, trAlign(ii));
    
    xbound(1) = min([xbound(1); XY(:,1)]);
    xbound(2) = max([xbound(2); XY(:,1)]);
    ybound(1) = min([ybound(1); XY(:,2)]);
    ybound(2) = max([ybound(2); XY(:,2)]);
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function im = getIm(name, tr)
MAXSIZE = 1024;

im = imread(name);

if size(im,3) > 1
    im = rgb2gray(im);
end
im = im2single(im);

szim = size(im);

if max(szim) > MAXSIZE
    im = imresize(im, MAXSIZE / max(szim));
end

if ~isempty(tr)
    im = applyTransformImage(im, tr);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function im = clipIm(im)
im(im > 1) = 1;
im(im < 0) = 0;
end
