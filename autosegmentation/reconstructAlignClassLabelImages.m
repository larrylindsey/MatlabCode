function reconstructAlignClassLabelImages(secdoc, cachedir,indices, ...
    classfun, labelfun)

if nargin < 5
    %labelfun = @(x, y) fijiTrainableSegmentationAnnotation(...
    %    im2double(x) > .75);\
    labelfun = @(x, y) better_segment(sizeImage(im2double(x), size(y)),...
        im2double(y));
    if nargin < 4
        classfun = @defaultClassFun;
        if nargin < 3
            cachedir = [pwd '/seg_cache'];
        end
    end
end

unix(['mkdir -p "' cachedir '"']);

[x y] = reconstructDomainBounds(secdoc);

secdoc = secdoc(indices);

parfor ii = 1:numel(secdoc)
    imclass = classfun(secdoc(ii));
    if ~isempty(imclass)        
        imclassTR = applyTransformImage(imclass, ...
            getTransform(secdoc(ii)), x, y);
        imwrite(imclassTR, sprintf('%s/class_%04d.png', cachedir, ...
            secdoc(ii).index));
        imclassTR = []; %#ok
        
        im = getIm(secdoc(ii));
        imTR = applyTransformImage(im, getTransform(secdoc(ii)), x, y);
        imwrite(imTR, sprintf('%s/image_%04d.png', cachedir, ...
            secdoc(ii).index));
        
        imlabel = uint16(labelfun(imclass, im));
        if ~isempty(imlabel)
            imlabelTR = applyTransformImage(imlabel, ...
                getTransform(secdoc(ii)), x, y, 'nearest');
            imwrite(imlabelTR, sprintf('%s/label_%04d.png', cachedir, ...
                secdoc(ii).index));
        end
    end
end

end

function tr = getTransform(secdoc)
imindex = secdoc(1).section.transImageIndex;
tr = secdoc(1).section.Transform(imindex);
end

function imclass = defaultClassFun(secdoc)
tr = getTransform(secdoc);
imfile = tr.Image.src;
imclass = imread(['class_' imfile]);
end

function im = getIm(secdoc)
tr = getTransform(secdoc);
imfile = tr.Image.src;
im = imread(imfile);
if size(im, 3) > 1
    im = rgb2gray(im);
end
end