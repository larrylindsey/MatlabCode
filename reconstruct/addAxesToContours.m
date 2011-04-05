function secdoc = addAxesToContours(secdoc, names, convex)

if nargin < 3
    convex = false;
end

transformTemplate = getTransTemplate(secdoc(1));

for is = 1:numel(secdoc)
    currSec = secdoc(is);
    newTrans = transformTemplate;
    for ic = 1:numel(names)
        
        contour = extractContour(currSec, names{ic});
        axisContours = repmat(struct, 0);
        for icc = 1:numel(contour)
            axisContours = cat(2, axisContours,...
                createAxisContours(contour(icc), convex));
        end
        
        if ic == 1
            newTrans.Contour = repmat(axisContours, 0);
        end
        newTrans.Contour = cat(2, newTrans.Contour, axisContours);
        secdoc(is).section.Transform(end + 1) = newTrans;
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function transformTemplate = getTransTemplate(sec)
transformTemplate = sec.section.Transform(1);
transformTemplate.dim = 0;
transformTemplate.xcoef = [0 1 0 0 0 0];
transformTemplate.ycoef = [0 0 1 0 0 0];
transformTemplate.T = [0 0; 1 0; 0 1; 0 0; 0 0; 0 0];
transformTemplate.Tinv = [0 0; 1 0; 0 1; 0 0; 0 0; 0 0];
transformTemplate.Contour = repmat(transformTemplate.Contour, 0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function axisContours = createAxisContours(contour, convex)
c = contourCentroid(contour);
c = c(1:2);
axisContours = rmfield(contour, {'section', 'z'});
axisContours = repmat(axisContours, [1 2]);

[spt d apt] = getMajorMinorAxes(c, contour.transPoints, convex);

axisContours(1).name = [contour.name 'AxisMajor'];
axisContours(1).closed = 0;
axisContours(1).border = [0 1 0];
axisContours(1).fill = [0 1 0];
axisContours(1).points = cat(1, spt(1,:), apt(1,:));
axisContours(1).transPoints = axisContours(1).points;

axisContours(2).name = [contour.name 'AxisMinor'];
axisContours(2).closed = 0;
axisContours(2).border = [1 0 0];
axisContours(2).fill = [1 0 0];
axisContours(2).points = cat(1, spt(2,:), apt(2,:));
axisContours(2).transPoints = axisContours(2).points;

end