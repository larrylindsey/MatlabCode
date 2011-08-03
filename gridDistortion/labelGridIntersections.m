function [rc_fit, model, model_err, L, crossEnergy] = ...
    labelGridIntersections(im, match, mask, crossEnergy, bwEnergy)
%Right angle triangle fit threshold
angle_thresh = .1 * pi;
%Ransac parameters
ransac_maxerr = .1;
ransac_inlier_factor = .25;
ransac_iter = 100;



if nargin < 3 || isempty(mask)
    mask = true(size(im));
end

if nargin <4
    tic;
    crossEnergy = normcorr(im, match);
    crossEnergy(not(mask)) = 0;
    clear imnormsqr;
    toc;
end

if nargin < 5
    
    tPercent = 0.99;
    
    [eh, ehx] = hist(crossEnergy(mask), 1000);
    pp = cumsum(eh);
    pp = pp / pp(end);
    thI = find(pp > tPercent, 1, 'first');
    p2 = pp(thI); x2 = ehx(thI);
    p1 = pp(thI - 1); x1 = ehx(thI - 1);
    
    tVal = ((tPercent - p1) / (p2 - p1)) * x2 + ...
        ((p2 - tPercent)/(p2 - p1)) * x1;

    clear crossEnergySorted;

    bwEnergy = crossEnergy > tVal;
end


L = sparse(bwlabel(bwEnergy, 4));
rc = getLabelPeaks(L);


[sample_space angles] = makeTriangleSampleSpace(rc);
sel = abs(angles - pi / 2) < angle_thresh;
sample_space = sample_space(sel,:);
ransac_min_inliers = round(max(L(:)) * ransac_inlier_factor);

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