function [cal, wcal] = calibrateIt(im, nomsize)

if ischar(im)
    im = imread(im);
end

if nargin < 2 || isempty(nomsize)
    nomsize = .463;
end


if isstruct(im)
    str = im;
else
    im = removeBlackBackground(im);
    
    str  = extractDistortionTransform(im);
end


% Get rc_found, scale from [-1 1] -> [0 1]
rc = str.rc_found + 1;
rc = rc / 2;

% Scale up to image coordinates
rc(:,1) = rc(:,1) * str.sz(1);
rc(:,2) = rc(:,2) * str.sz(2);

[vcal hcal] = doCalibration(rc, str.rc_grid);

figure;
plot(sort(nomsize ./ vcal), linspace(0, 1, numel(vcal)), 'r', ...
    sort(nomsize ./ hcal), linspace(0, 1, numel(hcal)), 'g', 'LineWidth', 2);
legend(sprintf('Vertical Calibration CDF, n = %d', numel(vcal)), ...
    sprintf('Horizontal Calibration CDF, n = %d', numel(hcal)));
xlabel('Units / pix');
ylabel('Fraction');
title('Calibration Statistics');

cal = nomsize / mean(cat(1, vcal, hcal));
n = numel(vcal) + numel(hcal);
wcal = nomsize / ((mean(vcal) * numel(hcal) + mean(hcal) * numel(vcal)) / n);

fprintf('Raw calibration: %0.4g units per pixel\n', cal);

fprintf('V/H norm calibration: %0.4g units per pixel\n', wcal);


% model(:,1) = model(:,1) * imsz(1) / 2;
% model(:,2) = model(:,2) * imsz(2) / 2;
% 
% cal = nomsize / mean(sqrt(sum(model.^2,1)));
% 
% fprintf('Calibration: %g units per pixel\n', cal);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vertcal, horizcal] = doCalibration(rc, grid)

nbd = cell(size(rc, 1), 1);

for ii = 1:size(rc, 1)
    [inbd isvert] = findNeighbors(grid, ii);
    localnbd = zeros(numel(inbd), 6);
    
    localnbd(:, 1:2) = grid(inbd, :);
    localnbd(:, 3:4) = repmat(grid(ii,:), [numel(inbd) 1]);
    localnbd(:, 6) = isvert;
    
    for jj = 1:numel(inbd)
        localnbd(jj, 5) = sum((rc(ii, :) - rc(inbd(jj), :)).^2);
    end
    nbd{ii} = localnbd;
end

nbd = cat(1, nbd{:});

vertsel = nbd(:,4) > 0;
horizsel = not(vertsel);

vertcal = sqrt(nbd(vertsel, 5));
horizcal = sqrt(nbd(horizsel, 5));

%cal(1:2,:) = sqrt(cal(1:2,:));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inbd isvert] = findNeighbors(grid, ii)
gi = grid(ii,:);

onbd = [-1 0; 0 -1; 0 1; 1 0];
% We're in row-column coordinates, so -1 0 is vertical, 0 -1 is horizontal.
isvert = [true; false; false; true];
inbd = zeros(size(onbd, 1), 1);


for jj = 1:size(onbd, 1)
    onj = onbd(jj,:);
    nb = gi + onj;
    sel = grid == repmat(nb, [size(grid, 1) 1]);
    sel = and(sel(:,1), sel(:,2));
    
    if any(sel)
        inbd(jj) = find(sel, 1, 'first');
        if inbd(jj) < ii
            inbd(jj) = 0;
        end
    end
end

isvert(inbd <= 0) = [];
inbd(inbd <= 0) = [];

end
