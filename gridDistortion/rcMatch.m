function [selA selB rc1 rc2] = rcMatch(rc1, rc2)

rc1 = setOrigin(rc1);
rc2 = setOrigin(rc2);

rc_offset = getOffset(rc1, rc2);
rc2 = rc2 - repmat(rc_offset, [size(rc2, 1), 1]);

[~, selA, selB] = intersect(rc1, rc2, 'rows');

end

function offset = getOffset(rc1, rc2)
im1 = gridData2im(rc1, ones(size(rc1,1),1));
im2 = gridData2im(rc2, ones(size(rc2,1),1));

im1(isnan(im1)) = 0;
im2(isnan(im2)) = 0;

im2o = round(size(im2) / 2);

oconv = imfilter(im1,im2);
mm = max(oconv(:));
[ox oy] = find(oconv == mm, 1, 'first');

offset = im2o - [ox oy];
end

function rc = setOrigin(rc)
rc = rc - repmat(min(rc, [], 1), [size(rc, 1) 1]);
end

