function im = montageImages(im1, tr1, im2, tr2)

bbox1 = [1 1; fliplr(size(im1))];
bbox2 = [1 1; fliplr(size(im2))];

if ~isempty(tr1)
    bbox1 = doTransform(bbox1, tr1);
else
    tr1 = ...
        regressionTransform([0 0; 1 0; 0 1; 1 1], [0 0; 1 0; 0 1; 1 1], 1);
    tr1 = rmfield(tr1, 'data');
end

if ~isempty(tr2)
    bbox2 = doTransform(bbox2, tr2);
else
    tr2 = ...
        regressionTransform([0 0; 1 0; 0 1; 1 1], [0 0; 1 0; 0 1; 1 1], 1);
    tr2 = rmfield(tr2, 'data');
end

boxMin = min(cat(1, bbox1, bbox2), [], 1);
boxMax = max(cat(1, bbox1, bbox2), [], 1);

bbox = cat(1, boxMin, boxMax);

im1tr = applyTransformImage(im1, tr1, bbox(:,1)', bbox(:,2)');
im2tr = applyTransformImage(im2, tr2, bbox(:,1)', bbox(:,2)');

mask1 = applyTransformImage(...
    true(size(im1)), tr1, bbox(:,1)', bbox(:,2)', 'nearest');
mask2 = applyTransformImage(...
    true(size(im2)), tr2, bbox(:,1)', bbox(:,2)', 'nearest');

norm = mask1 + mask2;
norm(norm == 0) = 1;

im = (im2double(im1tr) + im2double(im2tr)) ./ norm;