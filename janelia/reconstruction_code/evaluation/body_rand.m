function [p, r, diff] = body_rand(gt, as, params)
% [P, R] = BODY_RAND(GT, AS [, PARAMS])
% Compute the Rand precision and recall using a constant number of points
% per segmentation body.
% INPUT:
%   - GT: the ground truth segmentation
%   - AS: the automatic segmentation
%   - PARAMS: struct of options. Fields:
%      - perform_watershed_gt: if the ground truth is a boundary map
%      instead of a segmentation, run watershed on it. (default: true)
%      - perform_watershed_as: if the automatic segmentation is a boundary
%      map instead of a segmentation, run watershed on it. (default: true)
%      - as_ws_threshold: if the automatic segmentation boundary map has
%      continuous boundary values, threshold them by this threshold first.
%      (default: 1)
%      - markers_per_body: how many points to select in each body.
%      (default: 10)
%      - min_body_size: disregard bodies smaller than this size. (default:
%      1000)
%      - verbose: print runtime information
%      - debug: print debugging information
    if ~exist('params', 'var'),
        params = struct();
    end
    if ~isfield(params, 'perform_watershed_gt'),
        params.perform_watershed_gt = true;
    end
    if ~isfield(params, 'perform_watershed_as'),
        params.perform_watershed_as = true;
    end
    if ~isfield(params, 'as_ws_threshold'),
        params.as_ws_threshold = 1;
    end
    if ~isfield(params, 'markers_per_body'),
        params.markers_per_body = 10;
    end
    if ~isfield(params, 'min_body_size'),
        params.min_body_size = 1000;
    end
    if ~isfield(params, 'diff'),
        params.diff = false;
    end
    if ~isfield(params, 'debug'),
        params.debug = false;
    end
    
    if params.perform_watershed_gt,
        gt = watershed(gt, 6);
    end
    if params.perform_watershed_as,
        as = watershed(uint8(as >= params.as_ws_threshold), 6);
    end
    
    body_labels = sort(unique(gt));
    body_labels = body_labels(2:end); % ignore the 0 label
    body_label_counts = body_labels;
    for l=1:numel(body_labels),
        body_label_counts(l) = nnz(gt == body_labels(l));
    end
    body_labels = body_labels(body_label_counts >= params.min_body_size);
    n = numel(body_labels);
    m = params.markers_per_body;
    neg = nchoosek(n, 2) * m^2;
    pos = nchoosek(m, 2) * n;
    
    if ~isfield(params, 'submatrix_multiplier'),
        params.submatrix_multiplier = 2;
    end
    if ~isfield(params, 'submatrix_size'),
        params.submatrix_size = params.submatrix_multiplier*n;
    end
    
    gtmarkers = zeros([n*m 1], class(body_labels));
    asmarkers = zeros([n*m 1], class(body_labels));
    markerids = zeros([n*m 1]);
    
    for i=1:n,
        ids = randsample(find((gt == body_labels(i)) .* (as ~= 0)), m);
        markerids(((i-1)*m+1):(i*m)) = ids;
        gtmarkers(((i-1)*m+1):(i*m)) = gt(ids);
        asmarkers(((i-1)*m+1):(i*m)) = as(ids);
    end
    
    if params.debug,
        fprintf(['number of labels: %i,\n' ...
                'number of markers per label: %i,\n' ...
                'matrix size: %i,\n' ...
                'number of positive pairs: %i,\n' ...
                'number of negative pairs: %i,\n'], ...
                n, m, n*m, pos, neg);
    end
    
    s = params.submatrix_size;
    r = 0;
    p = 0;
    diff = struct('false_splits', zeros([s^2 2]), ...
                  'false_merges', zeros([s^2 2]));
    fs_idx = 1;
    fm_idx = 1;
    for i=1:(n*m/s),
        for j=i:(n*m/s),
            gtij = repmat(gtmarkers(((i-1)*s+1):(i*s)), [1 s]);
            gttij = repmat(gtmarkers(((j-1)*s+1):(j*s)), [1 s])';
            dgt = uint8(logical(gtij - gttij));
            asij = repmat(asmarkers(((i-1)*s+1):(i*s)), [1 s]);
            astij = repmat(asmarkers(((j-1)*s+1):(j*s)), [1 s])';
            das = uint8(logical(asij - astij));
            if params.diff,
                [rx, ry] = find((1-das).*dgt);
                [px, py] = find(das.*(1-dgt));
            end
            if i == j,
                r = r + nnz(das.*dgt)/2;
                p = p + (nnz((1-das).*(1-dgt)) - s)/2;
                if params.diff,
                    validp = px < py;
                    px = px(validp);
                    py = py(validp);
                    validr = rx < ry;
                    rx = rx(validr);
                    ry = ry(validr);
                end
            else
                r = r + nnz(das.*dgt);
                p = p + nnz((1-das).*(1-dgt));
            end
            if params.diff,
                nsplits = numel(px);
                nmerges = numel(rx);
                px = px + (i-1)*s;
                py = py + (j-1)*s;
                rx = rx + (i-1)*s;
                ry = ry + (j-1)*s;
                if fs_idx+nsplits-1 > size(diff.false_splits, 1),
                    diff.false_splits = vertcat(diff.false_splits, ...
                                        zeros(size(diff.false_splits)));
                end
                diff.false_splits(fs_idx:(fs_idx+nsplits-1),:) = ...
                    [asmarkers(px) asmarkers(py)];
                fs_idx = fs_idx + nsplits;
                if fm_idx+nmerges-1 > size(diff.false_merges, 1),
                    diff.false_merges = vertcat(diff.false_merges, ...
                                        zeros(size(diff.false_merges)));
                end
                diff.false_merges(fm_idx:(fm_idx+nmerges-1),:) = ...
                    [gtmarkers(rx) gtmarkers(ry)];
                fm_idx = fm_idx + nmerges;
            end
        end
    end
    diff.false_splits = uniquify(diff.false_splits);
    diff.false_merges = uniquify(diff.false_merges);
    r = r/neg;
    p = p/pos;
    
end

function m2 = uniquify(m)
    maxm = max(m(:));
    d = 10^(ceil(log10(maxm+1)));
    m = unique(m(:,1)+m(:,2)/d);
    if m(1) == 0,
        m = m(2:end);
    end
    m2 = zeros([numel(m) 2]);
    m2(:,1) = floor(m);
    m2(:,2) = rem(m, 1)*d;
    m2 = round(m2); % resolves approximation errors from floating point arithmetic
end