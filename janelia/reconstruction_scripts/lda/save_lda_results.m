% auc = zeros(0, 8);
% for i=1:8,
%     auc(i) = allperf{i}{4};
% end

datadir = '/groups/chklovskii/home/nuneziglesiasj/Projects/em_denoising/data/10x10x10_cropped/multi_image_tif/';
numsamps = [12500 25000 50000 100000];
for i=1:4,
    filename = sprintf('lda11_%i.tif', numsamps(i));
    imout = imfilter(double(im), allperf{i}{1}, 'symmetric');
    write_image_stack(imout, struct('filename', [datadir filename], 'histeq', 0, 'invert', 1));
    [testpred, testlabels] = sample_patches(imout(:,:,1:115), bmb(:,:,1:115), 50000, [1 1 1], struct('balanced', true));
    [rec, prec, ~, auc] = perfcurve(testlabels, testpred, 1, 'xCrit', 'reca', 'yCrit', 'prec');
    allperf{i}{2} = prec;
    allperf{i}{3} = rec;
    allperf{i}{4} = auc;
end