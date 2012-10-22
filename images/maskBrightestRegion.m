function mask = maskBrightestRegion(im, th)
%


% nevermind
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


mask = im >= (max(im(:)) - th);

% 
% if nargin < 3
%     if strcmp(class(im), 'uint8')
%         uth = max(im(:));
%         step = 1;        
%     else
%         uth = max(im(:));
%         step = (uth - lth) / 10;
%     end
% end
% 
% checkVal = lth:step:uth;
% 
% n = numel(checkVal);
% 
% masksum = zeros(size(im));
% 
% for ii = 1:n
%     mask = im >= checkVal(ii);
%     masksum = masksum + mask;
% end
% 
% mask = masksum > (n / 2);

end
