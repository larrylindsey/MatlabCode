function performance = lda_test(im, labels, params)
    if ~exist('params', 'var'),
        params = struct();
    end
    if ~isfield(params, 'patchsize'),
        params.patchsize = [11 11 11];
    end
    if ~isfield(params, 'samplesize'),
        params.samplesize = 150000;
    end
    if ~isfield(params, 'test_samplesize'),
        params.test_samplesize = params.samplesize;
    end
    if ~isfield(params, 'testz'),
        params.testz = 1:115;
    end
    if ~isfield(params, 'trainz'),
        params.trainz = 116:size(im, 3);
    end
    if ~isfield(params, 'labels_to_sample'),
        params.labels_to_sample = [0 1];
    end
    if ~isfield(params, 'zero_labels'),
        params.zero_labels = [];
    end
    if ~isfield(params, 'filename'),
        params.filename = '~/temp/lda.tif';
    end
    if ~isfield(params, 'histeq'),
        params.histeq = false;
    end
    im = double(im)/255;
    [patches, labelv] = sample_patches(im(:,:,params.trainz), ...
                            labels(:,:,params.trainz), params.samplesize, ...
                            params.patchsize, struct('balanced', true, ...
                            'label_list', params.labels_to_sample));
    zero_labels = params.zero_labels;
    if ~isempty(zero_labels),
        labels_to_zero = uint8(labelv == zero_labels(1));
        for i=2:numel(zero_labels),
            labels_to_zero = labels_to_zero .* (labelv == zero_labels(i));
        end
        labelv(logical(labels_to_zero)) = 0;
    end
    f = lda(patches, labelv, params.patchsize);
    p = patches * f(:);
    if (mean(p(labelv == 1)) - mean(p(labelv == 0))) < 0,
        f = -f;
    end
    imout = imfilter(im, f, 'symmetric');
    write_image_stack(imout, struct('filename', params.filename, ...
                                    'histeq', params.histeq, ...
                                    'invert', true));
    test_labels = labels;
    if ~isempty(zero_labels),
        labels_to_zero = uint8(test_labels == zero_labels(1));
        for i=2:numel(zero_labels),
            labels_to_zero = labels_to_zero .* (test_labels == zero_labels(i));
        end
        test_labels(logical(labels_to_zero)) = 0;
    end
    [testpred, testlabels] = sample_patches(imout(:,:,params.testz), ...
                                            test_labels, ...
                                            params.test_samplesize, ...
                                            [1 1 1], struct('balanced', true));
    [rec, prec, ~, auc] = perfcurve(testlabels, testpred, 1, ...
                                    'xCrit', 'reca', 'yCrit', 'prec');
    performance = {f, prec, rec, auc};
end