allperf = cell(1,8);
datadir = '/groups/chklovskii/home/nuneziglesiasj/Projects/em_denoising/data/10x10x10_cropped/multi_image_tif/';

numsamps = [12500 25000 50000 100000];
params = struct('labels_to_sample', [0 1], 'zero_labels', []);
fn_base = [datadir 'lda11_'];
for i=1:4,
    params.samplesize = numsamps(i);
    params.filename = [fn_base sprintf('%i', numsamps(i)) '.tif'];
    allperf{i} = lda_test(im, bmb, params);
end

numsamps = [18750 37500 75000 150000];
params = struct('labels_to_sample', [0 1 2], 'zero_labels', 2);
fn_base = [datadir 'lda11_mito_'];
for i=1:4,
    params.samplesize = numsamps(i);
    params.filename = [fn_base sprintf('%i', numsamps(i)) '.tif'];
    allperf{4+i} = lda_test(im, gtb, params);
end