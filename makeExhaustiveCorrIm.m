function [corrIm distIm] = makeExhaustiveCorrIm(im1, im2, res, ws, offset)

sz1 = size(im1);
sz2 = size(im2);

sz1 = sz1(1:2);
sz2 = sz2(1:2);

if nargin < 5
    if sz1 ~= sz2
        fprintf('No offset specified. Calculating one automagically\n');
    end    
    offset = round((sz1 - sz2) / 2);
    
    if nargin < 4

        if nargin < 3
            res = round(min([sz1 sz2] / 256));
        end
        ws = res * 4 + 1;        
    end
end

if ws / 2 == round(ws / 2)
    ws = ws + 1;
end

if isinteger(im1)
    im1 = mean(double(im1) / 255, 3);
end

if isinteger(im2)
    im2 = mean(double(im2) / 255, 3);
end

rowLim(1) = min([1 1 + offset(1)]);
rowLim(2) = max([sz1(1) 1 + offset(1)]);

colLim(1) = min([1 1 + offset(2)]);
colLim(2) = max([sz1(2) 1 + offset(2)]);

rowDom = rowLim(1):rowLim(2);
colDom = colLim(1):colLim(2);

nr = floor(numel(rowDom) / res);
nc = floor(numel(colDom) / res);

corrIm = zeros(nr, nc);
distIm = corrIm;

h = waitbar(0, 'Doot Doot');
iw = 0;
iwmax = nr * nc;

halfws = round(ws / 2);

%figure;
%is(1) = subplot(1, 2, 1);
%axis equal; colormap gray;
%is(2) = subplot(1, 2, 2);
%axis equal; colormap gray;
warning('off', 'MATLAB:divideByZero');
for i_r = 1:nr
    for i_c = 1:nc
        coord = ([i_r i_c] - 1) * res + 1;
        win1 = getWindow(im1, ws, coord, rowDom, colDom);
        win2 = getWindow(im2, ws * 3, coord - offset, rowDom, colDom);
%        imagesc(win1, 'Parent', is(1));
%        imagesc(win2, 'Parent', is(2));
%        keyboard;
        wincorr = imfilter(win1, win2);        
        maxCorr = max(wincorr(:));
        corrIm(i_r, i_c) = maxCorr;
        [locr locc] = find(wincorr == maxCorr);
        distvect = cat(2, locr, locc);
        dd = dist2(distvect, [halfws  halfws]);
        distIm(i_r, i_c) = min(dd);
        waitbar(iw / iwmax, h);
        iw = iw + 1;
    end
end
warning('on', 'MATLAB:divideByZero');

close(h);

corrIm(logical(isnan(corrIm))) = 0;

end

function win = getWindow(im, ws, coord, rd, cd)
win = zeros(ws);
sz = size(im);

winDomRows = rd(coord(1)) + (1:ws) - floor(ws / 2);
winDomCols = cd(coord(2)) + (1:ws) - floor(ws / 2);

rsel = and(logical(winDomRows > 0), logical(winDomRows <= sz(1)));
csel = and(logical(winDomCols > 0), logical(winDomCols <= sz(2)));

win(rsel, csel) = im(winDomRows(rsel), winDomCols(csel));
end


