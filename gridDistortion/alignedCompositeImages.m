function alignedCompositeImages(cache, imfiles, outdir)

unix(['mkdir -p ' outdir]);

% im1 = im2double(rgb2gray(imread(imfiles{1})));
ncpu = max(matlabpool('size'), 1);
splitImfiles0 = cell(1, ncpu);
splitImfiles1 = cell(1, ncpu);
splitFeat1 = cell(1, ncpu);
splitFeat2 = cell(1, ncpu);
splitIndices = cell(1, ncpu);

for ii = 1:ncpu
    splitIndices{ii} = ii:ncpu:numel(cache.feat1);
    splitFeat1{ii} = cache.feat1(splitIndices{ii});
    splitFeat2{ii} = cache.feat2(splitIndices{ii});
    splitImfiles0{ii} = imfiles(splitIndices{ii});
    splitImfiles1{ii} = imfiles(splitIndices{ii} + 1);
end


parfor cpu = 1:ncpu
    feat1 = splitFeat1{cpu};
    feat2 = splitFeat2{cpu};
    imfiles0 = splitImfiles0{cpu};
    imfiles1 = splitImfiles1{cpu};
    indices = splitIndices{cpu};
    
    for ii = 1:numel(feat1)
        im1 = im2double(rgb2gray(imread(imfiles0{ii})));
        im2 = im2double(rgb2gray(imread(imfiles1{ii})));
        tr = fitTransform(feat1{ii}, feat2{ii}, 1, @taylorMat);
        tr = fitInverseTransform(tr);
        tr.data = struct('n', 32, 'u', [-1 1], 'v', [-1 1]);
        im2tr = applyTransformImage(im2, tr, [-1 1], [-1 1]);
        %     imwrite(compositeImage(im1, im2tr), ...
        %         sprintf('%s/%03d aligned %03d.png', outdir, ii, ii + 1));
        imwrite(imresize(compositeImage(im1, im2tr), .25), ...
            sprintf('%s/th %03d aligned %03d.png', outdir, ...
            indices(ii), indices(ii) + 1));
%         fprintf('Just finished for section %d\n', ii);
%         im1 = im2;
    end
end
end
