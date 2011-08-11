function fh = ...
    pset1_part4(im , scale2, func, varargin)

if nargin < 3
    func = @hThenWResize;
end

%im = imread(filename);
im_dbl = double(im) / 255;

sz = size(im_dbl);
newsz = sz;
newsz(1:2) = round(newsz(1:2) .* scale2);

reduction = sz - newsz;

[im_seamed seam_path]= func(im, reduction(1:2), varargin{:});

seam_sz = size(im_seamed);

im_scaled = imresize(im_dbl, newsz(1:2));

scaled_sz = size(im_scaled);

fh = figure;

subplot(2, 2, 1);
image(im_dbl);
title(sprintf('Original Image, %g x %g', sz(2), sz(1)));

subplot(2, 2, 2);
image(im_seamed);
title(sprintf('Seam-Carved Image, %g x %g', seam_sz(2), seam_sz(1)));

subplot(2, 2, 3);
image(im_scaled);
title(sprintf('Conventionally Scaled Image, %g x %g', scaled_sz(2), ...
    scaled_sz(1)));

subplot(2,2,4);
ph = plot(seam_path(2,:), seam_path(1,:), 'k');
axis([0 sz(2) 0 sz(1)]);
title('Size History of Seam-Carved Image');
set(ylabel('Vertical Size'), 'FontSize', 14);
set(xlabel('Horizontal Size'), 'FontSize', 14);
set(ph, 'LineWidth', 2);

