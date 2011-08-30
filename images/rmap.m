function map = rmap(im1, im2, n)

if size(im1) ~= size(im2)
    error('Images must be the same size');
end

samp = 1:n;

map = zeros(size(im1) - n);

for rr = 1:size(map,1)
    for cc = 1:size(map,2)
        patch1 = im1((rr - 1) + samp, (cc - 1) + samp);
        patch2 = im2((rr - 1) + samp, (cc - 1) + samp);
        rmat = corrcoef(patch1(:), patch2(:));
        map(rr,cc) = rmat(2);
    end
end
