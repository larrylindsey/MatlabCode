function do_hw_1

seals = imread('seals.jpg');

groceries = imread('groceries.jpg');

trees = imread('trees.jpg');

skyline1 = imread('skyline1.jpg');

skyline2 = imread('skyline2.jpg');

raspberry = imread('raspberry.jpg');

trackstand = imread('trackstand.jpg');

imlist = {seals, groceries, trees, skyline1, skyline2, raspberry,...
    trackstand};

scalelist = [1 .8; .8 1; .5 .9; .9 .5];


% fh = do_part1(seals);
% drawnow;
% fh = cat(2, fh, do_part2(seals));
% drawnow;
% fh = cat(2, fh, do_part3(seals));
% drawnow;
% fh = cat(2,fh,do_part4(seals, 'Seals', groceries, 'Groceries', trees, 'Trees'));
% drawnow;
% 
% for i_f = 1:numel(fh)    
%     setTitle(fh(i_f));
%     file = sprintf('Fig_%g.jpg', i_f);
%     saveas(fh(i_f), file);
% end
% 
% %do_part3(seals);
% %do_part4(seals, 'Seals', groceries, 'Groceries');
% 
% for i_im = 1:numel(imlist)
%     for i_scale = 1:size(scalelist, 1)
%         fh = do_part5(imlist{i_im}, scalelist(i_scale,:),...
%             @greedyResize, @sumOfDerivativesEnergy);
%         setTitle(fh);
%         newname = sprintf('Fig5_%g_%g.jpg', i_im, i_scale);
%         saveas(fh, newname);
%         %close;
%     end
% end

for i_im = [4 5]
    for i_scale = [3 4]
        fh = do_part5(imlist{i_im}, scalelist(i_scale,:),...
            @hThenWResize, @sumOfDerivativesEnergy);
        setTitle(fh);
        newname = sprintf('Fig5_%g_%g_hw.jpg', i_im, i_scale);
        saveas(fh, newname);
        
        fh = do_part5(imlist{i_im}, scalelist(i_scale,:),...
            @wThenHResize, @sumOfDerivativesEnergy);
        setTitle(fh);
        newname = sprintf('Fig5_%g_%g_wh.jpg', i_im, i_scale);
        saveas(fh, newname);
    end

end

end

function fh = do_part1(seals)

seals_small = reduceWidth(seals, 325);

seals_small_letterbox = letterbox(seals_small, size(seals));

fh = figure;

subplot(1, 2, 1);
imshow(seals);
axis equal;
title('Original Seals Image');

subplot(1, 2, 2);
imshow(seals_small_letterbox);
axis equal;
title(sprintf('Seals Image With 325\nVertical Seams Removed'));

end

function fh = do_part2(seals)

fh = zeros(1, 2);

seals_energy = defaultEnergy(seals);

columnSeamCostmap = calculateCostmap(seals_energy);

rowSeamCostmap = calculateCostmap(seals_energy')';

fh(1) = figure;

subplot(2, 2, 1);
imshow(seals);
title('Seals Image');

subplot(2, 2, 2);
imagesc(seals_energy);
axis equal;
colorbar;
title('Seals Image Energy Function');

subplot(2, 2, 3);
imagesc(columnSeamCostmap);
axis equal;
colorbar;
title('Cost Map for Column Seams');

subplot(2, 2, 4);
imagesc(rowSeamCostmap);
axis equal;
colorbar;
title('Cost Map for Row Seams');

fh(2) = figure;


%Normalize for composite images
seals_energy = seals_energy / max(seals_energy(:));

cscEnergy = defaultEnergy(columnSeamCostmap);
cscEnergy = cscEnergy / max(cscEnergy(:));

rscEnergy = defaultEnergy(rowSeamCostmap);
rscEnergy = rscEnergy / max(rscEnergy(:));

im_composite_col = zeros(size(seals));
im_composite_col(:,:,1) = cscEnergy;
im_composite_col(:,:,2) = (cscEnergy + seals_energy) / 2;
im_composite_col(:,:,3) = seals_energy;

im_composite_row = zeros(size(seals));
im_composite_row(:,:,1) = rscEnergy;
im_composite_row(:,:,2) = (rscEnergy + seals_energy) /2;
im_composite_row(:,:,3) = seals_energy;


subplot(1, 2, 1);
imshow(im_composite_col);
title(sprintf('Composite Image\nof Column Cost Map Differential\nand Energy Function'));

subplot(1, 2, 2);
imshow(im_composite_row);
title(sprintf('Composite Image\nof Row Cost Map Differential\nand Energy Function'));
end

function fh = do_part3(seals)

fh = figure;

ee = defaultEnergy(seals);
cmcol = calculateCostmap(ee);
cmrow = calculateCostmap(ee')';

seamcol = calculate_seam(cmcol) + 1;
seamrow = calculate_seam(cmrow) + 1;

%Normalize for jet colormap
cmcol = cmcol * 64 / max(cmcol(:));
cmrow = cmrow * 64 / max(cmrow(:));

subplot(2, 2, 1);
imshow(seals);
ph(1) = display_seam(seamrow);
title(sprintf('Seals Image\nwith First Row Seam Indicated'));

subplot(2, 2, 2);
imshow(cmrow, jet);
drawnow;
ph(2) = display_seam(seamrow);
title(sprintf('Seals Image Row Seam Costmap\nwith First Row Seam Indicated'));

subplot(2, 2, 3);

imshow(seals);
ph(3) = display_seam(seamcol');
title(sprintf('Seals Image\nwith First Column Seam Indicated'));

subplot(2, 2, 4);
imshow(cmcol,jet);
drawnow;
ph(4) = display_seam(seamcol');
title(sprintf('Seals Image Row Seam Costmap\nwith First Row Seam Indicated'));

set(ph, 'LineWidth', 2);
set(ph, 'Color', [1 1 1]);

end

function fh = do_part4(varargin)

imCell = {varargin{1:2:end}};
imNames = {varargin{2:2:end}};

n = numel(imCell);

fh = zeros(1, n);

for i_im = 1:n
   im = imCell{i_im};
    
    sz = size(im);
    rc = sz(1:2);
    
    reduction = floor(rc(2) * .25);

    fh(i_im) = figure;
    
    subplot(2, 2, 1);
    imshow(im);
    title(sprintf(imNames{i_im}));
    
    imout_e1 = reduceWidth(im, reduction, @sumOfDerivativesEnergy);
    imout_e2 = reduceWidth(im, reduction, @derivativesOfGaussianEnergy);
    
    imout_grey_e1 = mean(imout_e1, 3);
    imout_grey_e2 = mean(imout_e2, 3);
    
    imout_compare = zeros(size(imout_e1));
    imout_compare(:,:,1) = imout_grey_e1;
    imout_compare(:,:,2) = (imout_grey_e1 + imout_grey_e2) / 2;
    imout_compare(:,:,3) = imout_grey_e2;
    
    subplot(2, 2, 2);
    imshow(letterbox(imout_compare, size(im)));
    title(sprintf(['Comparison of Reduced Images\nOrange Intensity due' ...
        ' to Left Image,\nCyan due to right']));
    
    subplot(2, 2, 3);
    imshow(letterbox(imout_e1, size(im)));
    title(sprintf('Image Reduced by\nSum of Magnitude of Derivatives'));
    
    subplot(2, 2, 4);
    imshow(letterbox(imout_e2, size(im)));
    title(sprintf(['Image Reduced by\nSum of Magnitude of Derivatives\n' ...
        'Smoothed by Gaussian h = 12, s = .5']));
    
end


end

function fh = do_part5(im, scale2, func, varargin)

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
imshow(im_dbl);
title(sprintf('Original Image, %g x %g', sz(2), sz(1)));

subplot(2, 2, 2);
imshow(letterbox(im_seamed, size(im_dbl)));
title(sprintf('Seam-Carved Image, %g x %g', seam_sz(2), seam_sz(1)));

subplot(2, 2, 3);
imshow(letterbox(im_scaled, size(im_dbl)));
title(sprintf('Conventionally Scaled Image,\n%g x %g', scaled_sz(2), ...
    scaled_sz(1)));

subplot(2,2,4);
ph = plot(seam_path(:,2), seam_path(:,1), 'k');
axis([0 sz(2) 0 sz(1)] + [-10 10 -10 10]);
grid on;
title('Size History of Seam-Carved Image');

set(ylabel('Vertical Size'), 'FontSize', 14);
set(xlabel('Horizontal Size'), 'FontSize', 14);
set(ph, 'LineWidth', 2);

end

function imout = letterbox(im, sz)

imout = zeros(sz);

rows = size(im, 1);
cols = size(im, 2);

dr = sz(1) - rows;
dc = sz(2) - cols;

ri = (1:rows) + floor(dr / 2);
ci = (1:cols) + floor(dc / 2);

imout(ri, ci, :) = im;

end

function setTitle(fh)
h_axes = get(fh, 'Children');
if iscell(h_axes)
    h_title = get(cat(1, h_axes{:}), 'Title');
else
    h_title = get(h_axes, 'Title');
end
set([h_title{:}], 'FontSize', 16);
end