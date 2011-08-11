function efunc = energy1(image_in_dbl)
image_in_gray = mean(image_in_dbl,3);
d_dx = [-1 1];
d_dy = [-1 ; 1];
efunc = zeros(size(image_in_gray));
efunc = efunc + abs(imfilter(image_in_gray, d_dx, 'same', 'replicate'));
efunc = efunc + abs(imfilter(image_in_gray, d_dy, 'same', 'replicate'));
