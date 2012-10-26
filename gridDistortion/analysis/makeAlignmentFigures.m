function [im1tr, im2tr, imAlign] = makeAlignmentFigures(id1, cache, imfiles)
id2 = id1 + 1;

trAlign = seriesAlignmentTransforms(cache, id1);
feat1 = cache.feat1{id1};
feat2 = cache.feat2{id1};
feat2tr = applyTransform(feat2, trAlign(id2));

im1 = im2single(imread(imfiles{id1}));
im2 = im2single(imread(imfiles{id2}));
im1 = mean(im1, 3);
im2 = mean(im2, 3);

XY = gridRC([-1 1], [-1 1]);
XYtr = applyTransform(XY, trAlign(id2));
x = cat(1, XY(:,1), XYtr(:,1));
y = cat(1, XY(:,2), XYtr(:,2));

x = [min(x), max(x)];
y = [min(y), max(y)];
px = [x(1) x(1) x(2) x(2) x(1)];
py = [y(1) y(2) y(2) y(1) y(1)];

if nargout > 0
    im1tr = applyTransformImage(im1, trAlign(id1), x, y);
    im2tr = applyTransformImage(im2, trAlign(id2), x, y);
    
    mask1tr = isnan(im1tr);
    mask2tr = isnan(im2tr);
    
    im1tr(mask1tr) = 0;
    im2tr(mask2tr) = 0;
    
    imAlign = (im1tr + im2tr) / 2;
    
    imAlign(and(mask1tr, mask2tr)) = 1;
    im1tr(mask1tr) = 1;
    im2tr(mask2tr) = 1;
end

imshow(im1, 'XData', [-1 1], 'YData', [-1 1]);
hold on;
plot(feat1(:,1), feat1(:,2), 'bx', 'LineWidth', 2);

figure;
imshow(im2, 'XData', [-1 1], 'YData', [-1 1]);
hold on;
plot(feat2(:,1), feat2(:,2), 'bo', 'LineWidth', 2);

figure;
plot(feat1(:,1), feat1(:,2), 'bx', feat2(:,1), feat2(:,2), 'bo', ...
    'LineWidth', 2);
hold on;
plot(cat(2, feat1(:,1), feat2(:,1))', cat(2, feat1(:,2), feat2(:,2))', 'k');
plot(px, py, 'w');
axis([x y]);
axis image
set(gca, 'XTickLabel', '', 'YTickLabel', '')

figure;
plot(feat1(:,1), feat1(:,2), 'bx', feat2tr(:,1), feat2tr(:,2), 'bo', ...
    'LineWidth', 2);
hold on;
plot(cat(2, feat1(:,1), feat2tr(:,1))', cat(2, feat1(:,2), feat2tr(:,2))', 'k');
plot(px, py, 'w');
axis([x y]);
axis image
set(gca, 'XTickLabel', '', 'YTickLabel', '')

end
