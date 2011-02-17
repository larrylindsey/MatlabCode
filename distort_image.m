function [im_out, mask] = distort_image(im, PP, coords)

x = linspace(-1,1,size(im,2));
y = linspace(-1,1,size(im,1));

[xx yy] = meshgrid(x,y);

pxx = zeros(size(xx));
pyy = zeros(size(yy));

for ii = 1:size(coords,2)
    pxx = pxx + PP(:,:,ii) * coords(1, ii);
    pyy = pyy + PP(:,:,ii) * coords(2, ii);
end

% im_out = zeros(size(im));
% 
% for ii = 1:size(im,3)
%     im_out(:,:,ii) = interp2(xx, yy, im(:,:,ii), pxx, pyy,'linear', 0);
% end

zz = ones(size(im));
for ii = 1:size(zz,3)
    zz(:,:,ii) = zz(:,:,ii) * ii;
end

rm_sz = [1 1 size(im,3)];

% xx = repmat(xx, rm_sz);
% yy = repmat(yy, rm_sz);
% pxx = repmat(pxx, rm_sz);
% pyy = repmat(pyy, rm_sz);
% 
% xx = single(xx);
% yy = single(yy);
% pxx = single(pxx);
% pyy = single(pyy);
% zz = single(zz);
% im = single(im);

%im_out = griddata3(pxx,pyy,zz,im,xx,yy,zz,'linear');
for ii = 1:3
    im_out(:,:,ii) = griddata_lite(pxx, pyy, im(:,:,ii), xx, yy);
end

im_out(im_out(:) > 1) = 1;
im_out(im_out(:) < 0) = 0;

% im_out = single(im_out);

mask = griddata(pxx(:,:,1), pyy(:,:,1), ones(size(im, 1), size(im, 2)),...
    xx(:,:,1), yy(:,:,1), 'nearest');
mask = repmat(mask, [1 1 size(im, 3)]);
% mask = single(mask);