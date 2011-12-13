function gridStruct = findRoughLines(im, rStr, spt)

showDisplay = true;

if nargin < 3
    spt = 0;
end

if spt < 1

    gaussFilt = fspecial('gaussian', 50, 9);
    edgeFiltVert = fspecial('sobel');
    edgeFiltHoriz = edgeFiltVert';


    edgeFiltVertBlur = imfilter(gaussFilt, edgeFiltVert);
    edgeFiltHorizBlur = imfilter(gaussFilt, edgeFiltHoriz);
    imVertBlur = imfilter(im, edgeFiltVertBlur);
    imHorizBlur = imfilter(im, edgeFiltHorizBlur);


    edgeMap = (abs(imHorizBlur) + abs(imVertBlur));
    
    %cropMask = autocropGrid(edgeMap);
    cropMask = true(size(edgeMap));
    
    %edgeMapSort = sort(edgeMap(cropMask));
    edgeMapSort = sort(edgeMap(:));
    edgeMapT = edgeMapSort(round(7/ 8 * numel(edgeMapSort)));

    bwEdge = and(logical(edgeMap > edgeMapT), cropMask);
    
%     bwEdge = bwareaopen(bwEdge, 1024, 4);
%     bwEdge = imerode(bwEdge, strel('disk', 4));

else
    disp('Skipping Edge Detection...');
    bwEdge = rStr.bw;
    cropMask = rStr.cropMask;
end

if spt < 2
    [h t r] = hough(bwEdge, 'ThetaResolution', .1);
else
    disp('Skipping Hough Transform...');
    h = rStr.h;
    t = rStr.t;
    r = rStr.r;
end    

% Estimate roughly, the number of lines we expect to find.
% Each grating cell should be about 128px on a side, give or take a factor
% 2 or so.
nEstLines = round(2 *max(size(im)) / 128);

h = medfilt2(h, [3 3]);
pp = findPeaks2(h);
hpeaks = h(pp);
hpeaks_sorted = sort(hpeaks);
sel = (h .* pp > hpeaks_sorted(round(end/2)));

pksel = and(sel, h > hpeaks_sorted(end - nEstLines));
[pr pc] = find(pksel);
peaks = cat(2, pr, pc);

angleDistribution = sum(sel .* h, 1) ./ sum(sel, 1);
angleDistribution = angleDistribution(1:ceil(end / 2)) + ...
    angleDistribution((floor(end/2) + 1):end);
tDistribution = t((floor(end / 2) + 1):end);

% % Experimentally determined normalization.
% % This is used to remove a periodic bias for peaks in the hough space with
% % regard to line orientation
% 
% fnorm = (13 - cosd(4 * tDistribution)) / 13;
% angleDistribution = angleDistribution ./ fnorm;

angleDistribution(isnan(angleDistribution)) = 0;
angleDistributionFT = fft(angleDistribution);
angleDistributionFT([2 end]) = 0;
angleDistribution = ifft(angleDistributionFT);

[~, iFortyFive] = min(abs(tDistribution - 45));
angleDistribution(iFortyFive) = 0;

[~, imax] = max(angleDistribution);
deg = tDistribution(imax);

%fprintf('Rotating Grid');

%gridSupport = imrotate(true(size(im)), deg, 'nearest', 'crop');
%gridOut = imrotate(im, deg, 'bilinear', 'crop');
%gridOut(not(gridSupport)) = median(im(:));

if showDisplay
    close all;
    setShowDisplay(figure);
    imshow(bwEdge);
end

fprintf('Done');

gridStruct.deg = deg;
gridStruct.bw = bwEdge;
gridStruct.peaks = peaks;
gridStruct.h = h;
gridStruct.t = t;
gridStruct.r = r;
gridStruct.sz = size(im);
gridStruct.cropMask = cropMask;
