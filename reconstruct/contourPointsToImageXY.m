function contours = contourPointsToImageXY(secdoc, contours)

if nargin < 2
    contours = enumerateContours(secdoc);
end

if isstruct(contours)
    contours = {contours};
end

if ischar(contours{1})
    for ii = 1:numel(contours)
        contours{ii} = extractContour(secdoc, contours{ii});
    end
end

secIndices = [secdoc.index];

for ii = 1:numel(contours)
    contour = contours{ii};
    
    for jj = 1:numel(contour)
        sec = secdoc(secIndices == contour(jj).section).section;
        tr = sec.Transform(sec.transImageIndex);
        
        tpts = contour(jj).transPoints;
        tpts = applyInverseTransform(tpts, tr);
        tpts = tpts / sec.Transform(sec.transImageIndex).Image.mag;
        imsz = max(sec.Transform(sec.transImageIndex).Contour(1).points, [], 1);
        imszblock = repmat(imsz, [size(tpts, 1), 1]);

        contour(jj).pixelPts = tpts;
        contour(jj).pixelPtsSc = 2 * tpts ./ imszblock - 1;
        contour(jj).imsz = imsz;        
    end
    contours{ii} = contour;
end

end
