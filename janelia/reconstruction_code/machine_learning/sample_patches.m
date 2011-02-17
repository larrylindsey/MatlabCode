function [patches, patch_labels] = sample_patches(im, labels, numsamples, ...
                                                         patchsize, params)
% [PATCHES, PATCH_LABELS] = SAMPLE_PATCHES(IM, LABELS, NUMSAMPLES,
% PATCHSIZE, PARAMS(BALANCED, LABEL_LIST, TEST_PERCENT))
    if any(mod(patchsize, 2)==0),
        error('sample_patches: patchsize should be odd in each dimension.');
    end
    if nargin < 5,
        params = struct();
    end
    if ~isfield(params, 'balanced'),
        params.balanced = true;
    end
    if ~isfield(params, 'label_list'),
        params.label_list = [0 1];
    end
    if ~isfield(params, 'test_percent'),
        params.test_percent = [0.5, 0.5];
    else
        if sum(params.test_percent(:)) ~= 1,
            error('params.test_percent should add up to 1.');
        end
    end
    
    % Calculate start and end positions (stay away from block faces)
    offset = floor(patchsize/2);
    minpos = offset+1;
    minx = minpos(1);
    miny = minpos(2);
    minz = minpos(3);
    maxpos = size(im)-offset;
    maxx = maxpos(1);
    maxy = maxpos(2);
    maxz = maxpos(3);
    
    startz = zeros([numel(params.test_percent), 1]);
    endz = zeros(size(startz));
    startz(1) = 1;
    endz(1) = params.test_percent(1)*size(labels,3);
    for i=2:numel(startz),
        startz(i) = endz(i-1)+1;
        endz(i) = endz(i-1)+params.test_percent(i)*size(labels,3);
    end
    if params.balanced,
        numsamples_per_label = floor(numsamples/numel(params.label_list));
        numsamples = numsamples_per_label*numel(params.label_list);
    end
    patches = zeros([numsamples, prod(patchsize)]);
    patch_labels = zeros([numsamples, 1]);
    for i=1:numel(params.label_list),
        cur_label = params.label_list(i);
        [x y] = find(labels == cur_label);
        z = ceil(y/size(im,2));
        y = mod(y, size(im,2));
        selector = find( ...
                    (x >= minx) .* (x <= maxx) .* ...
                    (y >= miny) .* (y <= maxy) .* ...
                    (z >= minz) .* (z <= maxz) ...
                    );
        population = [x(selector), y(selector), z(selector)];
        n = size(population, 1);
        k = min([n numsamples_per_label]);
        sample = randsample(n, k);
        ids = population(sample, :);
        for j=1:numsamples_per_label,
            k = (i-1)*numsamples_per_label+j;
            top = ids(j,:)-offset;
            bot = ids(j,:)+offset;
            patch = im(top(1):bot(1), top(2):bot(2), top(3):bot(3));
            patches(k,:) = patch(:);
            patch_labels(k) = labels(ids(j,1), ids(j,2), ids(j,3));
        end
        
    end
end