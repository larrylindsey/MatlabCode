function secdoc = addContoursTo(secdoc, contours, indices, template)
% function secdoc = addContoursTo(secdoc, contours, indices, template)
%
% Adds contours to the reconstruct series represented in secdoc
%
% Inputs
% secdoc - struct representation of a Reconstruct series, as from
%   readReconstruct
% contours - {n}, cell array containing contours stored in [p 2], where p
%   is the number of points.
% indices - [m] where m <= n. Typically m should equal n. For each k in 
%   [1..m], indices(k) indicates the section in which contours{k} should be
%   placed.
% template - may be one of:
%   [struct], in which case template must contain the following fields:
%       name - the contour's name
%       hidden - [logical], indicates whether the new contour will be
%          hidden in Reconstruct.
%       closed - [logical], indicates whether the contour is a closed or
%          open contour.
%       simplified - [logical], just set it to true. Whatevs.
%       border [1 3] - rgb triplet in [0..1], indicates border color.
%       fill [1 3] - rgb triplet in [0..1], indicates fill color.
%       mode [1] - I can't remember what this does. Seems like 11 is a good
%           choice.
%   {[char array]}, in which case the first cell should contain the name of
%       an existing contour to use as a template, and the second should 
%       contain the new name.


if numel(contours) ~= numel(indices)
    warning('addContoursTo:cardinality', ...
        'contours and indices contain differing numbers of elements');
end

if iscell(template)
    sections = [secdoc.section];
    transforms = [sections.Transform];
    it = 1;
    dontHaveIt = true;
    name = template{1};
    while it < numel(transforms) && dontHaveIt
        transContours = transforms(it).Contour;
        for ic = 1:numel(transContours)
            if strcmp(name, transContours(ic).name)
                templStruct = transContours(ic);
                dontHaveIt = false;
            end
        end
        it = it + 1;
    end
    
    templStruct.name = template{2};
    template = templStruct;
end

[x y] = reconstructDomainBounds(secdoc);

transformTemplate = makeTransformTemplate(x,y);

secIndex = [secdoc.index];

for ii = 1:numel(indices)
    secSel = find(secIndex == indices(ii), 1, 'first');
    if isempty(secSel)
        warning('addContoursTo:sectionNotFound',...
            'No section found with index %d, skipping', indices(ii));
    else
        transform = transformTemplate;
        transform.Contour = template;
        transform.Contour.points = contours{ii};
        transform.Contour.transPoints = contours{ii};
        secdoc(secSel).section.Transform = catstruct(...
            secdoc(secSel).section.Transform, transform);
    end
end
end

function transformTemplate = makeTransformTemplate(x, y)
transformTemplate.dim = 0;
transformTemplate.xcoef = [0 1 0 0 0 0];
transformTemplate.ycoef = [0 0 1 0 0 0];
transformTemplate.Contour = [];
transformTemplate.Image = [];
transformTemplate.type = @taylorMat;
transformTemplate.T = zeros(6, 2);
transformTemplate.T([2 9]) = 1;
transformTemplate.Tinv = transformTemplate.T;
transformTemplate.order = 2;
transformTemplate.iDim = 2;
transformTemplate.oDim = 2;
transformTemplate.data.n = 32;
transformTemplate.data.u = x;
transformTemplate.data.v = y;
end