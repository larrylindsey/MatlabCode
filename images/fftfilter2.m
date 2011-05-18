function out = fftfilter2(im, kern)

if size(im, 3) > 1 
    error('This function supports only two dimensions of cross corr.');
end

imrc = size(im);

krc = size(kern);
krc = krc([1 2]);

lbrc = imrc + krc;

if lbrc(1) == 2 * round(lbrc(1) / 2)
    lbrc(1) = lbrc(1) + 1;
end
if lbrc(2) == 2 * round(lbrc(2) / 2)
    lbrc(2) = lbrc(2) + 1;
end


select_out = logical(letterbox(true(size(im)), lbrc));

imlb = letterbox(im, lbrc);
kernlb = letterbox(kern, lbrc);

im_ft = fft2(imlb);
%kern_ft = fft2(flipud(fliplr(kernlb)));
nk = size(kern, 3);

outlb = zeros([lbrc nk]);

for i_k = 1:nk
    kern_ft = fft2(kernlb(:,:,i_k));
    out_ft = im_ft.*conj(kern_ft);
    outlb(:,:,i_k) = ifft2(out_ft);
end

%Shift outlb around so that we have matches in the right place
rshift = mod((1:lbrc(1))  + round(lbrc(1) / 2), lbrc(1));
rshift(logical(rshift == 0)) = lbrc(1);

cshift = mod((1:lbrc(2))  + round(lbrc(2) / 2), lbrc(2));
cshift(logical(cshift == 0)) = lbrc(2);


out = repmat(im(1), [size(im) nk]);

for i_k = 1:nk
    out_temp = outlb(rshift, cshift, i_k);
    out(:, :, i_k) = reshape(out_temp(select_out), size(im));
end

if isreal(im) && isreal(kern)
    out = real(out);
end


end