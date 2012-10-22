function mask = maskBrightestRegion(im, th)
%

if size(im, 3) > 1
    im = mean(im, 3);
end


mask = im >= (max(im(:)) - th);

end
