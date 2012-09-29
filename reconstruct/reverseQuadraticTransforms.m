function secdoc = reverseQuadraticTransforms(secdoc, doID)

if nargin < 2
    doID = false;
end

for ii = 1:numel(secdoc)
    sec = secdoc(ii);
    imIndex = sec.section.transImageIndex;
    imtr = sec.section.Transform(imIndex);
    
    if imtr.dim > 3                
        % Linearize imtr
        ltr = imtr;
        if doID
            ltr.dim = 0;
            ltr.xcoef = [0 1 0 0 0 0];
            ltr.ycoef = [0 0 1 0 0 0];
            ltr.T = [0 0; 0 1; 1 0];
            ltr.order = 1;
            ltr.monomialOrder = [0 0; 0 1; 1 0];
            ltr = fitInverseTransform(ltr);
        else
            ltr.dim = 3;
            ltr.xcoef(4:end) = 0;
            ltr.ycoef(4:end) = 0;
            ltr.T(4:6,:) = 0;
            ltr = fitInverseTransform(ltr);
            ltr.xcoef(1:3) = ltr.inv.T([1 3 2], 1);
            ltr.ycoef(1:3) = ltr.inv.T([1 3 2], 2);
        end
        
        % Iterate through all the Transforms in this section
        for it = 1:numel(sec.section.Transform)
            % Don't mess with the image transform subfields
            if it ~= imIndex                
                trans = sec.section.Transform(it);
                
                % Iterate through all the Contours in this Transform
                for ic = 1:numel(trans.Contour)
                    contour = trans.Contour(ic);
                    
                    % We deal only with pre-transformed points
                    pts = contour.transPoints;
                    pts = applyInverseTransform(pts, imtr);
                    pts = applyTransform(pts, ltr);
                    
                    % Now that they're "fixed", put them back.
                    contour.points = pts;
                    contour.transPoints = pts;
                    
                    trans.Contour(ic) = contour;
                end
                
                % Set trans to identity
                trans.xcoef = [0 1 0 0 0 0];
                trans.ycoef = [0 0 1 0 0 0];
                trans.T = [0 0; 0 1; 1 0];
                trans.order = 1;
                trans.monomialOrder = [0 0; 0 1; 1 0];
                trans = fitInverseTransform(trans);
                
                % write back
                sec.section.Transform(it) = trans;
            end
        end
        
        sec.section.Transform(imIndex) = ltr;
        secdoc(ii) = sec;
    end
end


end
