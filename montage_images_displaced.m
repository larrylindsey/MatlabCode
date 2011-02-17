function [imlt confidence xshift yshift] = ...
    montage_images_displaced(imlt, imrt)

imlt = single(imlt);
imrt = single(imrt);

[xshift yshift] = edge_align_images(imlt, imrt);

[hlt wlt] = size(imlt);
[hrt wrt] = size(imrt);

% Calculate final montage dimensions

montage_w = wlt + wrt - xshift -1;
hmin = min(hlt, hrt);
hmax = max(hlt, hrt);


if yshift > 0
    montage_h = max(hlt, hrt) + yshift - 1;
else
    h_over = max(hmax - hmin - yshift,0);
    montage_h = hmax + h_over - 1;
end

%imlt_merge = zeros(montage_h, montage_w);
%imrt_merge = imlt_merge;

% Calculate locations of old images within new montage images

% Left image
%leftx = 1:wlt;
lefty_beg = - min(0, yshift);
%lefty = lefty_beg + (1:hlt);

lth = size(imlt ,1);
ltwadd = montage_w - size(imlt,2);
imlt = cat(2, imlt, single(zeros(lth, ltwadd)));

ltw = size(imlt, 2);
lth = size(imlt, 1);
lthadd1 = lefty_beg;
lthadd2 = montage_h - lefty_beg - lth + 1;
imlt = cat(1, single(zeros(lthadd1, ltw)),...
    imlt, single(zeros(lthadd2, ltw)));

%imlt_merge(lefty, leftx) = imlt;

%Right image
%rightx = (1:wrt) + wlt - xshift - 1;
righty_beg = max(yshift, 0);
%righty = righty_beg + (1:hrt);
%imrt_merge(righty, rightx) = imrt;

rth = size(imrt, 1);
rtwadd = montage_w - size(imrt, 2);
imrt = cat(2, single(zeros(rth, rtwadd)), imrt);

rtw = size(imrt, 2);
rth = size(imrt, 1);
rthadd1 = righty_beg;
rthadd2 = montage_h - righty_beg + 1 - rth;
imrt = cat(1, single(zeros(rthadd1, rtw)), imrt,...
    single(zeros(rthadd2, rtw)));

%keyboard
%Make the merge rectangle
lt_control = zeros(1,4);
lt_control(1) = 1;

if yshift == 0
    lt_control(4) = -1;
elseif yshift > 0
    lt_control(4) = 1;
else
    lt_control(4) = 0;
end

lt_control(3) = 0;

if hrt + yshift - hlt == 0
    lt_control(2) = -1;
elseif hrt + yshift - hlt < 0
    lt_control(2) = 1;
else
    lt_control(2) = 0;
end

xmin = wlt - xshift;
xmax = wlt;
    
ymin = abs(yshift) + 1;
ymax = ymin + hmin - 1;

lt_merge_rect = blend_rectangle([xmin ymin], [xmax ymax], lt_control);
rt_merge_rect = 1 - lt_merge_rect;

xdom = xmin:xmax;
ydom = ymin:ymax;

confidence = corr2(imlt(ydom,xdom), imrt(ydom,xdom));

imlt(ydom, xdom) = imlt(ydom, xdom).*lt_merge_rect;
imrt(ydom, xdom) = imrt(ydom, xdom).*rt_merge_rect;

imlt = imlt + imrt;
