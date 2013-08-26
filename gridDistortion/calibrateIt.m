function cal = calibrateIt(im, nomsize, imsz)

if ischar(im)
    im = imread(im);
end

if nargin < 2 || isempty(nomsize)
    nomsize = .463;
end


if numel(im) == 4
    if nargin < 3
        error('Need image size to go with model');
    end
    model = im;
else
    im = removeBlackBackground(im);
    
    imsz = size(im);
    
    [~, ~, ~, model] = extractDistortionTransform(im);
end

model(:,1) = model(:,1) * imsz(1) / 2;
model(:,2) = model(:,2) * imsz(2) / 2;

cal = nomsize / mean(sqrt(sum(model.^2,1)));

fprintf('Calibration: %g units per pixel\n', cal);

end
