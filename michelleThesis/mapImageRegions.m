function mapImageRegions(rootFolder)

global c_k c_mito s_inf s_sup s_mito
c_k = [255  0	0;
    255	255	0;
    255  0    243;
    0    0   0] / 255;

c_mito = [0 0 1];
%infimum
s_inf = [2500, 50000, 250000];
% supremum
s_sup = [s_inf(2), s_inf(3), inf];
s_mito = 256;

if nargin < 1
    rootFolder = './';
end

if rootFolder(end) ~= '/' && rootFolder(end) ~= '\'
    rootFolder = [rootFolder '/'];
end

% List the files in rootDir, then select only the ones that are directories
dd = dir(rootFolder);
dirsel = [dd.isdir] == 1;
dd = dd(dirsel);
extA = 'noMit.tif';
extM = 'wMit.tif';

if rootFolder(end) ~= '/' && rootFolder(end) ~= '\'
    rootFolder = [rootFolder '/'];
end

% List the files in rootDir, then select only the ones that are directories
dd = dir(rootFolder);
dirsel = [dd.isdir] == 1;
dd = dd(dirsel);

allAFiles = {};
allMFiles = {};

for ii = 1:numel(dd)
    if ~(strcmp(dd(ii).name, '.') || strcmp(dd(ii).name, '..'))
        dirname = [rootFolder dd(ii).name];
        adir = dir([dirname '/*' extA]);
        mdir = dir([dirname '/*' extM]);
        
        if numel(adir) ~= numel(mdir)
            error(['Found %d annotation files and %d mito files for' ...
                ' dir %s.'], numel(adir), numel(mdir), dirname);
        else
            
            afiles = strcat([dirname '/'], sort({adir.name}));
            mfiles = strcat([dirname '/'], sort({mdir.name}));
            
            allAFiles = cat(2, allAFiles, afiles);
            allMFiles = cat(2, allMFiles, mfiles);
        end
    end
end


for i_f = 1:numel(allAFiles)
    mapImageHelper(allAFiles{i_f}, allMFiles{i_f});
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mapImageHelper(aFile, mFile)
global c_k s_inf s_sup c_mito s_mito
im_a = im2single(imread(aFile));
im_m = im2single(imread(mFile));
im_a = im_a(:,:,1:3);
im_m = im_m(:,:,1:3);
z = zeros(size(im_a, 1), size(im_a, 2), 'single');
w_a = z;
mask_m = false(size(z));
k = size(c_k, 1);

w = [.33 .67 1];

for i_k = 1:k
    mask = getMask(im_a, c_k(i_k,:));
    cc = bwconncomp(mask);
    for i_n = 1:cc.NumObjects
        px = cc.PixelIdxList{i_n};
        n_px = numel(px);
        if n_px > s_inf(1)
            if i_k < 3
                i_w = find(s_inf < n_px & s_sup >= n_px);
                w_a(px) = w(i_w(1));
            else
                w_a(px) = 1;
            end
        end
    end
end

mito_mask = getMask(im_m, c_mito);
cc = bwconncomp(mito_mask);
for i_n = 1:cc.NumObjects
    px = cc.PixelIdxList{i_n};
    if numel(px) > s_mito
        mask_m(px) = true;
    end
end

im_a2 = im_a .* repmat(w_a, [1 1 3]);
im_m2 = im_a2;
im_m2(repmat(mask_m, [1 1 3])) = im_m(repmat(mask_m, [1 1 3]));

imwrite(im2uint8(im_a2), newName(aFile));
imwrite(im2uint8(im_m2), newName(mFile));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mapName = newName(name)
[a, b, c] = fileparts(name);
if isempty(a)
    prefix = b;
else
    prefix = [a '/' b];
end

mapName = [prefix '.map' c];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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