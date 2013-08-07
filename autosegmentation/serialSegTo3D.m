function bwl = serialSegTo3D(bwl)
% bwl = serialSegTo3D(bwl)
%   Converts a serial 2D segmentation into a 3D segmentation.

bwl = double(bwl);

n = size(bwl, 3);

m = max(max(bwl(:,:,1)));

for ii = 2:n
    slice = bwl(:,:,ii);
    
    zmap = slice == 0;
    slice = slice + m;
    slice(zmap) = 0;
    
    bwl(:,:,ii) = slice;
    
    m = max(slice(:));
end
