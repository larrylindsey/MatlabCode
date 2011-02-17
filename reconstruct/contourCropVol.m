function contourCropVol(name, secdoc, scale, aval)

if nargin < 3
    scale = 1;
end

if nargin < 4
    aval = .5;
end

contour = extractContour(secdoc, name);

dz = secdoc(1).section.thickness;

csec = [contour.section];

usec = sort(unique(csec));

figure; hold on; colormap gray;

for i_sec = 1:numel(usec)
    currSec = usec(i_sec);
    
    secIndex = find([secdoc.index] == currSec, 1);

    [im x y] = getImage(secdoc, secIndex, scale);
    [X Y] = meshgrid(x, y);

    z = currSec * dz;
    Z = ones(size(X)) * z;

    contourIndex = find(csec == currSec);

    for i_c = 1:numel(contourIndex)
        [contourIm box r c] = contourSelect(im, contour(i_c), X, Y);
        sh = surface(reshape(X(box), [r c]),...
            reshape(Y(box), [r c]),...
            reshape(Z(box), [r c]), contourIm, 'LineStyle', 'none');
        alpha(sh, aval);
    end
    drawnow; if i_sec == 1; axis equal; view(45, 45); end
end

end

function [contourIm box r c] = contourSelect(im, contour, X, Y)

inContour = inpolygon(X, Y, contour.transPoints(:,1),...
    contour.transPoints(:,2));

csel = find(logical(sum(inContour, 1) > 0));
rsel = find(logical(sum(inContour, 2) > 0));

r = numel(rsel);
c = numel(csel);

box = false(size(X));
box(rsel, csel) = true;

inContour = reshape(inContour(box), [r c]);
contourIm = im(rsel, csel) .* inContour;
contourIm(~inContour) = nan;

end

function [imtr x y] = getImage(recstr, secindex, scale)

section = recstr(secindex).section;
imTransStr = section.Transform(section.transImageIndex);
im = imTransStr.Image;
imcontour = imTransStr.Contour;

fwfunc = doCubicTransform2(false);
revfunc = doCubicTransform2(true);


Ttot = cat(1, imTransStr.T, imTransStr.Tinv);
tr = maketform('custom', 2, 2, fwfunc, revfunc, Ttot);

imdata = im2double(rgb2gray(imread(im.src)));

if scale ~= 1
    imdata = imresize(imdata, scale, 'bilinear');
end

udata = [0 max(imcontour.points(:,1)) * im.mag];
vdata = [0 max(imcontour.points(:,2)) * im.mag];
[imtr x y]= imtransform(flipud(imdata), tr, 'UData', udata, 'VData', vdata);

x = linspace(x(1), x(2), size(imtr, 2));
y = linspace(y(1), y(2), size(imtr, 1));

end