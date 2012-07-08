function contours = applyTransformContourPixels(contours, tr, secdoc, force)
% contours = applyTransformContourPixels(contours, tr, [force])
%
%  Applies an image-level transform to a contour struct.
%
% contours - cell array of struct array
% tr - transform, as returned by fitTransform
% force - boolean, set to true to force a clean calculation. In other
%         words, forces the calculation of pixelPts, etc.

if nargin < 4
    force = false;
end

if force || ~isfield(contours{1}, 'pixelPts')
    contours = contourPointsToPixels(secdoc, contours);
end

secIndices = [secdoc.index];

for ii = 1:numel(contours)
    contour = contours{ii};
    for jj = 1:numel(contour)
        sec = secdoc(secIndices == contour(jj).section).section;
        recTr = sec.Transform(sec.transImageIndex);                
        
        tpts = contour(jj).pixelPts;
        
        n = size(tpts, 1);        
        imszblock = repmat(contour(jj).imsz, [n 1]);
        
        tpts = 2 * tpts ./ imszblock - 1;
        tpts = applyTransform(tpts, tr);
        tpts = imszblock .* (tpts + 1) / 2;        
        contour(jj).pixelPtsTr = tpts;
        tpts = tpts * sec.Transform(sec.transImageIndex).Image.mag;
        
        contour(jj).transPointsTr = applyTransform(tpts, recTr);
        
    end   
    contours{ii} = contour;
end

end
