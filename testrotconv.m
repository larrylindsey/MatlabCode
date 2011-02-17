function [corr imcorr] = testrotconv(d)
sz = 1000;
noise = randn(sz);
angles = 0:23;
corr = zeros(size(angles));

cropbeg = floor(sz / 4);
cropend = ceil(3 * sz / 4);
cropstep = ceil((cropend - cropbeg) / 500);
d = d * cropstep;
noise = conv2(noise, ones(d)/d/d,'same');
noisecrop = noise(cropbeg:cropstep:cropend,cropbeg:cropstep:cropend);
noisecropft = fft2(noisecrop);
imcorr = zeros([size(noisecrop) length(angles)]);

for i_angle = 1:length(angles)
    fprintf('%d\n',i_angle);
    angle = angles(i_angle);
    noiserot = imrotate(noise, angle, 'nearest', 'crop');
    noiserot = conv2(noiserot, ones(d)/d/d,'same');
    noiserot = noiserot(cropbeg:cropstep:cropend, cropbeg:cropstep:cropend);
    noiserotft = fft2(noiserot);
    noiseconv = ifft2(noiserotft .* conj(noisecropft));
    corr(i_angle) = max(max(noiseconv([1:2 end-(0:1),1:2 end-(0:1)])));
    imcorr(:,:,i_angle) = noiseconv;
end



end