function cal = calibrateIt(im, nomsize)

if ischar(im)
    im = imread(im);
end

if nargin < 2
    nomsize = .463;
end

[~, ~, ~, model] = extractDistortionTransform(im);

cal = nomsize / (mean(sqrt(model(:,1).^2 + model(:,2).^2)) * size(im,1)/2);

fprintf('Calibration: %g units per pixel\n', cal);

end
