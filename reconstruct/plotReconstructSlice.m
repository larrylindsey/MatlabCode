function fig = plotReconstructSlice(recstr, index, scale)

%fwfunc = doCubicTransform2(false);
%revfunc = doCubicTransform2(true);

secindex = find([recstr.index] == index, 1);

if isempty(secindex)
    error('Could not find requested index, %d', index);
end

section = recstr(secindex).section;
imTransStr = section.Transform(section.transImageIndex);
im = imTransStr.Image;
imcontour = imTransStr.Contour;

%Ttot = cat(1, imTransStr.T, imTransStr.Tinv);



imdata = imread(im.src);
if size(imdata, 3) > 1
    imdata = rgb2gray(imdata);
end

imdata = im2double(imdata);

if nargin > 2
    imdata = imresize(imdata, scale, 'bilinear');
end

if ~isfield(imTransStr, 'data') || ~isfield(imTransStr.data, 'u') ...
        || ~isfield(imTransStr.data, 'v');
    imTransStr.data.u = [0 max(imcontour.points(:,1)) * im.mag];
    imTransStr.data.v = [0 max(imcontour.points(:,2)) * im.mag];
end


[imtr x y] = applyTransformImage(imdata, imTransStr);
imtr = flipud(imtr);

%imtransform(flipud(imdata), tr, 'UData', udata, 'VData', vdata);

x = linspace(x(1), x(2), size(imtr, 2));
y = linspace(y(1), y(2), size(imtr, 1));


fig = figure;

set(fig, 'UserData', index);

imagesc(x, y, imtr); axis xy; colormap gray; axis image;
hold on;

for i_t = 1:numel(section.Transform)
    if i_t ~= section.transImageIndex
        contours = section.Transform(i_t).Contour;
        for i_c = 1:numel(contours)
            pts = contours(i_c).transPoints;
            if contours(i_c).closed
                pts = cat(1, pts, pts(1,:));
            end
            plot(pts(:,1), pts(:,2), 'LineWidth', 2, 'Color',...
                contours(i_c).border);
        end
    end
end
