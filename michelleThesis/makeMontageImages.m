function maxDim = makeMontageImages(rootFolder, animalIdAndGroup,...
    bestMontage, maxDim)

gid_a = animalIdAndGroup(:,2);
aid_a = animalIdAndGroup(:,1);

[~, filenames] = getMichelleMapImages(rootFolder);

gid_f = zeros(size(filenames)) - 1;


% There's probably a matlab-y way to do this that i've missed
for ii = 1:numel(aid_a)
    aidstr = num2str(aid_a(ii));
    locs = strfind(filenames, aidstr);
    for i_f = 1:numel(locs)
        if ~isempty(locs{i_f})
            gid_f(i_f) = gid_a(ii);
        end
    end
    
end


% group id in 0..5
youngsel = gid_f < 3;
oldsel = not(youngsel);

ncol = 10;

if nargin < 4
    maxDim = getMaxDim(filenames);
end

bestFilenames = getBest(filenames, bestMontage);

imBest = makeMontageImage(bestFilenames, 6, 1, 300, maxDim);
imOld = makeMontageImage(filenames(oldsel), ncol, 1, 300, maxDim);
imYoung = makeMontageImage(filenames(youngsel), ncol, 1, 300, maxDim);
imwrite(imBest, 'Best-montage.png');
imwrite(imOld, 'Old-montage.png');
imwrite(imYoung, 'Young-montage.png');

end

function imMontage = makeMontageImage(fileNames, ncol, marginInches, dpi, maxDim)


borderWeight = 2;

totWidthPx = (8.5 - marginInches) * dpi + borderWeight;
imSidePx = floor(totWidthPx / ncol) - borderWeight;
nrow = ceil(numel(fileNames) / ncol);
montageRC = [nrow ncol] * imSidePx;
imMontage = zeros([montageRC 3]);


scaleFactor = imSidePx / maxDim;

figure;
subplot(1, 2, 1);
mh = imshow(imMontage);
subplot(1, 2, 2);
iah = gca;
ih = imshow(zeros(imSidePx));

for ii = 1:numel(fileNames)
    file = fileNames{ii};
    if ~isempty(file)
        im = im2double(imread(file));
        im = im(:,:,1:3);
        im = imresize(im, scaleFactor);
        
        im = rotateCapillary(im);
        
        im(im > 1) = 1;
        im(im < 0) = 0;
        
        set(ih, 'CData', im);
        title(iah, fileNames{ii});
        
        C = mod(ii, ncol);
        if C == 0
            C = ncol;
        end
        R = ceil(ii/ ncol);
        
        c = (1:size(im,2)) + (C - 1) * imSidePx;
        r = (1:size(im,1)) + (R - 1) * imSidePx;
        imMontage(r, c, :) = im;
        set(mh, 'CData', imMontage);
        drawnow;
    end
end
end


function imOut = rotateCapillary(im)


e = -2;
imOut = [];

for i_d = 1:4
    imr = imrotate(im, (i_d - 1) * 90);
    sz = size(imr);
    [R, C] = meshgrid(linspace(0, 1, sz(2)), linspace(0, 1, sz(1)));
    okEnergy  = 1.2 * R + .8 * C;
    
    
    [ctr_r, ctr_c] = getCapCenter(imr);
    eNew = okEnergy(ctr_r, ctr_c);
    if eNew > e
        imOut = imr;
        e = eNew;
    end
end

end

function [r, c] = getCapCenter(im)
cap_color = .75 * [255  0    243] / 255;
mask = getMask(im, cap_color);
[r, c] = find(mask);
r = round(mean(r));
c = round(mean(c));
end

function mask = getMask(im, c)
% Creates a mask over im for a given color

% Start out everywhere true.
mask = true(size(im,1), size(im,2));
for ic = 1:numel(c)
    % Set the mask to false if the image's channel
    % doesn't match the corresponding component of
    % the color. This works for RGB or grayscale.
    mask = and(mask, abs(im(:,:,ic) - c(ic)) < (16 / 256));
end


end

function d = getMaxDim(fileNames)

d = 0;

for i_f = 1:numel(fileNames)
    sz = size(imread(fileNames{i_f}));
    dlast = d;
    d = max([d sz]);
    if d > dlast
        fprintf('%s was bigger than %d\n', fileNames{i_f}, dlast);
    end
end


end

function bestfiles = getBest(filenames, bestMontage)
b = size(bestMontage, 1);
f = numel(filenames);
bestfiles = {1, b};

for i_b = 1:b
    a_id = bestMontage{i_b, 2};
    exp = bestMontage{i_b, 1};
    for i_f = 1:f
        filename = filenames{i_f};
        if ~isempty(strfind(filename, a_id))
            if ~isempty(strfind(filename, exp))
                bestfiles{i_b} = filename;
            end
        end
    end
end

end