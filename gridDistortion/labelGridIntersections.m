function [rc_fit, model, model_err, L, crossEnergy] = ...
    labelGridIntersections(im, match, mask, crossEnergy, bwEnergy)
%Right angle triangle fit threshold
angle_thresh = .1 * pi;
%Ransac parameters
ransac_maxerr = .1;
ransac_inlier_factor = .25;
ransac_iter = 100;
expectAngle = pi / 2;



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
rc = getLabelPeaks(L, crossEnergy);

rc = simpleNonMaxSuppression(rc, crossEnergy, sqrt(2) * mean(size(match))/2);


[sample_space angles] = makeTriangleSampleSpace(rc, expectAngle);
% sel = abs(angles - expectAngle) < angle_thresh;
% sample_space = sample_space(sel,:);
ransac_min_inliers = round(size(rc,1) * ransac_inlier_factor);

[fit_space, model, model_err] = ransac(sample_space,...
    @createTriangleModel, @measureTriangleModel, ransac_maxerr,...
    ransac_min_inliers, 1, ransac_iter, rc);
    
rc_fit = rc(fit_space(:,5), :);

if getShowDisplay
    f = figure;
    imagesc(crossEnergy); axis image;
    hold on;
    ph(1) = plot(rc(:,2), rc(:,1), 'bx', 'LineWidth', 2);
    ph(2) = plot(rc_fit(:,2), rc_fit(:,1), 'g+', 'LineWidth', 2);
    set(f, 'UserData', ph);
end

% model in indicator vectors model(1,:), model(2,:)
v1 = model(1,:);
v2 = model(2,:);
% Find indicator that best corresponds to e1
a1 = atan2(v1(2), v1(1)) / pi;
a2 = atan2(v2(2), v2(1)) / pi;
a1m = mod(a1, 1);
a2m = mod(a2, 1);
if a1m > .25 && a1m <= .75
    if a2m <= .25 || a2m > .75
        model = model([2 1], :);
        %vtemp = v1;
        v1 = v2;
        %v2 = vtemp;
    end
end

if v1(1) < 0
    model(1,:) = -model(1,:);
end

if det(model) < 0
    model(2,:) = -model(2,:);
end

end

function [rc lh] = getLabelPeaks(L, crossEnergy)

n = full(max(L(:)));
rc = zeros(n, 2);
lh = zeros(n, 1);
r = cell(n, 1);
c = r;
ceSamp = r;

for ii = 1:n
    [r{ii} c{ii}] = find(L == ii);
    ceSamp{ii} = crossEnergy(sub2ind(size(crossEnergy), r{ii}, c{ii}));
end

parfor ii = 1:n
    [~, imax] = max(ceSamp{ii}); 
    rc(ii,:) = [r{ii}(imax), c{ii}(imax)];
    %rc(ii,:) = [mean(r) mean(c)];
    lh(ii) = numel(r{ii});
end

end

