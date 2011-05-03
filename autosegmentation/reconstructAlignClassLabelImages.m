function reconstructAlignClassLabelImages(secdoc, cachedir, ...
    classfun, labelfun)

if nargin < 4
    labelfun = @(x) fijiTrainableSegmentationAnnotation(...
        im2double(x) > .75);
    if nargin < 3
        classfun = @defaultClassFun;
        if nargin < 2
            cachedir = [pwd '/seg_cache'];
        end
    end
end

unix(['mkdir -p "' cachedir '"']);

[x y] = reconstructDomainBounds(secdoc);

parfor ii = 1:numel(secdoc)
    imclass = classfun(secdoc(ii));
    if ~isempty(imclass)        
        imclassTR = applyTransformImage(imclass, ...
            getTransform(secdoc(ii)), x, y);
        imwrite(imclassTR, sprintf('%s/class_%04d.png', cachedir, ...
            secdoc(ii).index));
        imclassTR = []; %#ok
        
        imlabel = uint16(labelfun(imclass));
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