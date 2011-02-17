function [cname dist] = selectContour(fig, secdoc)

if isempty(find(get(0, 'Children') == fig, 1))
    index = fig;
    fig = plotReconstructSlice(secdoc, fig, .25);
else
    index = get(fig, 'UserData');
end

figure(fig);
sec = secdoc(index).section;

title('Click inside contour to select');
[x y] = ginput(1);
title('');
dist = [];
cname = {};
for i_t = 1:numel(sec.Transform)
    if i_t ~= sec.transImageIndex
        currTran = sec.Transform(i_t);
        for i_c = 1:numel(currTran.Contour)
            currContour = currTran.Contour(i_c);
            if inpolygon(x, y, currContour.transPoints(:,1), ...
                    currContour.transPoints(:,2))
                d = rms(currContour.transPoints(1,:) - [x y]);
                cname = {cname{:} currContour.name};
                dist = [dist d];
            end
        end
    end
end
