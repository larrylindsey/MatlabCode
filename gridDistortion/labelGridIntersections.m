function [rc_fit, model, model_err, L, crossEnergy] = ...
    labelGridIntersections(im, match, mask, crossEnergy, bwEnergy)
%

%function [bwLabel crossEnergy rowLabel columnLabel match rowModel colModel] = ...
%    labelGridIntersections(im, match, crossEnergy)

%Initialize data stuffs
%strelSz = 6;
%se = strel('disk', strelSz);

%Right angle triangle fit threshold
angle_thresh = .1 * pi;
%Ransac parameters
ransac_maxerr = .1;
ransac_inlier_factor = .25;
ransac_iter = 100;


%[rr cc] = meshgrid(1:sz(1), 1:sz(2));
% 
% if nargin < 2 || (nargin <3 && isempty(match))
%     sz = size(im);
%     match = round(mean(sz(1:2)) / 20);
% end

% if numel(match) == 1 && nargin < 3
%     match = round(match / 2);
%     imshow(im);
%     title('Click on a good grid crossing, please');
%     [r c] = ginput(1);
%     r = round(r);
%     c = round(c);
%     close;
%     match = im(c + (-match:match), r + (-match:match));
% end


if nargin < 3 || isempty(mask)
    mask = true(size(im));
end

if nargin <4
    tic;

    %imnormsqr = fftfilter2(im.*im, ones(size(match)));
    %crossEnergy = fftfilter2(im, match) ./ sqrt(imnormsqr);

%     if false%max(size(im)) > 10240
%         %This takes EONS, but uses less memory, so will actually
%         %complete
%         fprintf('Creating Correlation Map\n');
%         imnormsqr = imfilter(im.*im, ones(size(match)));
%         fprintf('Finished norm map\n');
%         crossEnergy = imfilter(im, im2double(match)) ./ sqrt(imnormsqr);
%         fprintf('Finished corr map\n');
%     else
        imnormsqr = fftfilter2(im.*im, single(ones(size(match))));
        imnormsqr(imnormsqr < 0) = 0;
        crossEnergy = fftfilter2(im, match) ./ sqrt(imnormsqr);
        crossEnergy = crossEnergy .* mask;
%     end

    clear imnormsqr;
    toc;
end

if nargin < 5
    crossEnergySorted = sort(crossEnergy(mask));
    tVal = crossEnergySorted(round(.99 * numel(crossEnergySorted)));

    clear crossEnergySorted;

    %bwEnergy = imopen(crossEnergy > tVal, se);
    bwEnergy = crossEnergy > tVal;
end

%imwrite(bwEnergy, 'bwEnergy.png');
% 
% rHist = sum(bwEnergy, 1);
% cHist = sum(bwEnergy, 2);
% 
% [columnLabel colModel] = getRegularPeaks(cHist);
% [rowLabel rowModel] = getRegularPeaks(rHist');
% rowLabel = rowLabel';
% 
% bwLabel = bwlabel(bwEnergy, 8);
% 
% bwLabel = bwLabel .* repmat(logical(columnLabel > 0), size(rowLabel)) .* ...
%     repmat(logical(rowLabel > 0), size(columnLabel));

L = sparse(bwlabel(bwEnergy, 4));
[rc lh] = getLabelPeaks(L);

lh_sigma = sqrt(var(lh));
sel = lh < (mean(lh) + lh_sigma);
rc = rc(sel,:);

[sample_space angles] = makeTriangleSampleSpace(rc);
sel = abs(angles - pi / 2) < angle_thresh;
sample_space = sample_space(sel,:);
ransac_min_inliers = round(max(L(:)) * ransac_inlier_factor);
%rc = rc(sample_space(:,5),:);

[fit_space, model, model_err] = ransac(sample_space,...
    @createTriangleModel, @measureTriangleModel, ransac_maxerr,...
    ransac_min_inliers, 1, ransac_iter, rc);
    
rc_fit = rc(fit_space(:,5), :);

end

function [rc lh] = getLabelPeaks(L)

n = full(max(L(:)));
rc = zeros(n, 2);
lh = zeros(n, 1);

parfor ii = 1:n
    [r c] = find(L == ii);
    rc(ii,:) = [mean(r) mean(c)];
    lh(ii) = numel(r);
end

end