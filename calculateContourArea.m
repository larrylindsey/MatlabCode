function contours = calculateContourArea(contours, trAlign)

for ii = 1:numel(trAlign)
    trAlign(ii).secIndex = trAlign(ii).secIndex(1);
end
trSecIndex = [trAlign.secIndex];

if numel(unique(trSecIndex)) ~= numel(trSecIndex)
    error('Found multiple transforms for some sections');
end

for ii = 1:numel(contours)
    % contours is a cell index of arrays of struct
    % so, contour is now a struct array
    contour = contours{ii};
    for jj = 1:numel(contour)
        secSel = trSecIndex == contour(jj).section;
        if any(secSel)
            tra = trAlign(secSel);
            pix = contour(jj).pixelPtsSc;
            
            if ~isempty(tra.preTr)
                pix = applyTransform(pix, tra.preTr);
            end
            
            pix = applyTransform(pix, tra);
            pix(:,1) = (pix(:,1) + 1) * contour(jj).imsz(1) / 2;
            pix(:,2) = (pix(:,2) + 1) * contour(jj).imsz(2) / 2;
            pix = pix * contour(jj).mag;
            
            contour(jj).area = polyarea(pix(:,1), pix(:,2));
            contour(jj).slabVol = contour(jj).area * contour(jj).thickness;
            contour(jj).areaCalc = true;
            contour(jj).transPointsTr = pix;
        else
            contour(jj).area = [];            
            contour(jj).areaCalc = false;
            contour(jj).slabVol = [];
            contour(jj).transPointsTr = [];
        end        
    end
    contours{ii} = contour;
end
