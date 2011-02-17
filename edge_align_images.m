function [xshift yshift] = edge_align_images(imlt, imrt, guessx, guessy)
ltsize = size(imlt);
rtsize = size(imrt);

smallsize = min(ltsize, rtsize);
xsm = smallsize(2);
ysm = smallsize(1);

if nargin < 4
    guessy = round(ysm * .25);
    if nargin < 3
        guessx = round(xsm / 4);
    end
end

guessx = abs(guessx);
guessy = abs(guessy);

vertselect = (1:guessy) - ceil(guessy/2);

%ltend = imlt(guessy:(end+1-guessy),end + 1 - fliplr(1:guessx));
%rtend = imrt(guessy:(end+1-guessy),1:guessx);

ltend = imlt(vertselect + floor(end/2), end + 1 - fliplr(1:guessx));
rtend = imrt(vertselect + floor(end/2), 1:guessx);

ltendft = fft2(ltend);
rtendft = fft2(rtend);

imcorr = ifft2(ltendft.*conj(rtendft));

[corry corrx] = size(imcorr);

[yshift xshift] = find(imcorr == max(imcorr(:)));

%xshift = mod(xshift + floor(corrx / 2), corrx) - floor(corrx / 2);
xshift = corrx - xshift;
yshift = mod(yshift - 1 + floor(corry / 2), corry) - floor(corry / 2);

yshift = yshift + floor((size(imlt,1) - size(imrt,1))/2);