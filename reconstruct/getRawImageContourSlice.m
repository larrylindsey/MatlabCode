function slice = getRawImageContourSlice(secdoc, index, cmatch, bwfun)

if anyControlChars(cmatch)
    doExact = false;
else
    doExact = true;
end

if nargin < 4
    bwfun = @identity;
end

[trans imcontour mag secindex] = getParms(secdoc, index);

xmax = max(imcontour.points(:, 1));
ymax = max(imcontour.points(:, 2));
slice = false(ymax, xmax);

x = single(1:xmax);
y = single(1:ymax);
y = y(end:-1:1);
[X Y] = meshgrid(x, y);

sec = secdoc(secindex);

cEnum = enumerateContours(sec);
if doExact
    contours = {cEnum{strmatch(cmatch, cEnum, 'exact')}};
else
%     re = regexp(cEnum, cmatch);
%     sel = false(1, numel(cEnum));
%     for i_r = 1:numel(re)
%         if ~isempty(re{i_r})
%             sel(i_r) = true;
%         end
%     end
%     contours = {cEnum{sel}};
    contours = strsearch(cEnum, cmatch);
end

h = waitbar(0, 'doing some stuff');
M = numel(contours);
m = 1 / M;
q = 0;
for i_c = 1:M
    thisContour = extractContour(sec, contours{i_c});
    
    N = numel(thisContour);
    n = 1 / N;
    for ii = 1:N
        slice = getSlice(slice, trans, mag, X, Y, thisContour(ii), bwfun);
        waitbar(q + m * (ii * n), h);
    end
    q = i_c * m;
    waitbar(q, h);
end
close(h);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trans imcontour mag secindex] = getParms(secdoc, index)

secindex = find([secdoc.index] == index, 1);

if isempty(secindex)
    error('Could not find requested index, %d', index);
end

section = secdoc(secindex).section;

trans = section.Transform(section.transImageIndex);
mag = trans.Image.mag;
imcontour = trans.Contour;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slice = getSlice(slice, trans, mag, X, Y, contour, bwfun)

pts = contour.transPoints;

pts = doTransform(pts, trans, true);

pts = pts ./ mag;

if contour.closed    
    slice = slice | bwfun(inpolygon(X, Y, pts(:,1), pts(:,2)));
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eh = anyControlChars(pat)
cchars = '[^*$|<>\';
eh = false;
ii = 1;
while ~eh && ii < numel(cchars)
    eh = any(pat == cchars(ii));
    ii = ii + 1;
end
end