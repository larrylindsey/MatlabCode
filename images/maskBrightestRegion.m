function mask = maskBrightestRegion(im, lth, uth, step)
% mask = maskBrightestRegion(im, lth, uth, step)
%   im - the image to mask
%   lth - the lower threshold to check
%   uth - the upper threshold to check
%   step - the step between lth
%
% Masks an image by checking against progressively larger thresholds.
% Returns the median mask.
%
%
% mask = median(blockmask, 3), with blockmask(:,:,ii) = im >= checkTH(ii)

if size(im, 3) > 1
    im = mean(im, 3);
end

if nargin < 3
    if strcmp(class(im), 'uint8')
        uth = max(im(:));
        step = 1;        
    else
        uth = max(im(:));
        step = (uth - lth) / 10;
    end
end

checkVal = lth:step:uth;

n = numel(checkVal);

mask = false([size(im) n]);

for ii = 1:n
    mask(:,:,ii) = im >= checkVal(ii);
end

mask = median(mask, 3);

end
