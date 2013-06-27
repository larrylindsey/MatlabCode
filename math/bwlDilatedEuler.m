function e = bwlDilatedEuler(bwl3d)

s = size(bwl3d,3);
n = max(bwl3d(:));
e = zeros(n, s);



parfor jj = 1:s
    slice = bwl3d(:,:,jj);
    
    ep = zeros(n,1);
    seq = unique(slice(:));
    seq(seq == 0) = [];
    
    for ii = 1:numel(seq)
        kk = seq(ii);
        mask = slice == kk;
        ep(kk) = dEulerNum(mask, 15);
    end
    
    e(:,jj) = ep;
end

e = sum(e,2);
end

function e = dEulerNum(mask, s)
mask = imdilate(mask, strel('disk', s));
rp = regionprops(mask, 'EulerNumber');
if isfield(rp, 'EulerNumber')
    e = sum([rp.EulerNumber]);
else
    e = 0;
end
end
