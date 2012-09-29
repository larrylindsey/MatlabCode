function secdoc = unityTransformContours(secdoc)

for ii = 1:numel(secdoc)
    sec = secdoc(ii);
    imIndex = sec.section.transImageIndex;
    
    % Iterate through all the Transforms in this section
    for it = 1:numel(sec.section.Transform)
        % Don't mess with the image transform subfields
        if it ~= imIndex
            trans = sec.section.Transform(it);
            
            % Iterate through all the Contours in this Transform
            for ic = 1:numel(trans.Contour)
                contour = trans.Contour(ic);
                contour.points = contour.transPoints;
                trans.Contour(ic) = contour;
            end
            
            % Set trans to identity
            trans.xcoef = [0 1 0 0 0 0];
            trans.ycoef = [0 0 1 0 0 0];
            trans.T = [0 0; 1 0; 1 0];
            trans.order = 1;
            trans.monomialOrder = [0 0; 1 0; 1 0];
            trans = fitInverseTransform(trans);
            
            % write back
            sec.section.Transform(it) = trans;
        end
    end
    
    secdoc(ii) = sec;
end
end

