function rc = removeSmallCC(rc, n)

imrc = gridData2im(rc, ones([size(rc, 1) 1]));
imrc(isnan(imrc)) = 0;

L = bwlabel(imrc > 0, 4);
nl = labelHistogram(L);

for ii = 1:numel(nl)
    if nl(ii) < n
        imrc(L == ii) = -1;
    end
end

keyboard
