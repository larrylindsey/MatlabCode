function e = bwlDilatedEuler(bwl3d)

n = max(bwl3d(:));
e = zeros(n, 1);



parfor ii = 1:n    
    mask = bwl3d == ii;
    [~,~,ij] = find(mask, 1, 'first');
    mask = mask(:,:,ij);
    e(ii) = dEulerNum(mask, 15);
end

end

function e = dEulerNum(mask, s)
mask = imdilate(mask, strel('disk', s));
rp = regionprops(mask, 'Euler');
e = rp.EulerNumber;
end
