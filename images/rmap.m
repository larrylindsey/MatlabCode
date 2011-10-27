function map = rmap(im1, im2, n)

ncpu = matlabpool('size');

if ncpu < 1
    ncpu = 1;
end

if size(im1) ~= size(im2)
    error('Images must be the same size');
end

if size(im1, 3) > 1
    im1 = rgb2gray(im1);
end

if size(im2, 3) > 1
    im2 = rgb2gray(im2);
end

im1 = im2single(im1);
im2 = im2single(im2);

if numel(n) > 1
    m = n(2);
    n = n(1);
else
    m = n;
end

samp = 1:n;

%outmap = cell(1, ncpu);

%parfor ii = 1:ncpu

%localmap = zeros(size(im1) - n);
rr = 1:m:(size(im1,1) - n);
cc = 1:m:(size(im1,2) - n);

map = zeros(numel(rr), numel(cc));

for irr = 1:numel(rr)
    for icc = 1:numel(cc)
        patch1 = im1((rr(irr) - 1) + samp, (cc(icc) - 1) + samp);
        patch2 = im2((rr(irr) - 1) + samp, (cc(icc) - 1) + samp);
        rmat = corrcoef(patch1(:), patch2(:));
        map(irr,icc) = rmat(2);
    end
end

%map{ii} = localmap;
%    outmap{ii} = localmap;
%end

% 
% map = zeros(size(im1) - n);
% for ii = 1:ncpu
%     map = map + outmap{ii};
% end