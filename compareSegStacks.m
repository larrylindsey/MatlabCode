function cmapout = compareSegStacks(bwl1, bwl2, cmap)
n = max(size(bwl1, 3), size(bwl2, 3));
k = max(max(bwl1(:)), max(bwl2(:)));


if nargin < 3
    cmap = hsv(k + 1);
    cmap = cmap(randperm(size(cmap, 1)),:);
    cmap(1,:) = 0;
end

for ii = 1:n
    doImage(1, bwl1, cmap, ii);
    doImage(2, bwl2, cmap, ii);
    pause;
end

if nargout > 0
    cmapout = cmap;
end

end

function doImage(s, bwl, cmap, ii)
subplot(1, 2, s);
imagesc(bwl(:,:,ii));
axis image;
colormap(cmap);
caxis([0 size(cmap,1)-1]);
title(sprintf('Slice %d', ii));
end
