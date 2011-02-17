function [pr_table, colnames] = all_bwr(spec)
% [PR_TABLE, COLNAMES] = ALL_BWR(SPEC)
% Generate a table of body-wise rand precision recall values encompassing 
% multiple segmentation methods, potentially from different parts of an
% image stack.
% INPUT:
%  SPEC: a structure with the following fields. 
%    - ws_filenames: cell array of the filenames for the watershed files
%    - amb_file_patterns: cell array of the amb files corresponding to the
%    watershed files
%    - ts: cell array of the thresholds for each amb_file_pattern OR single
%    set of thresholds, identical for each amb_file_pattern
%    - gt: ground truth segmentation
%    - index_sets: which parts of the image to evaluate
%    - colname_base: the base name of each column
%    - trim: the amount of margin to remove from evaluation in BWR
%    - min_body_size: the smallest body (in voxels) to be considered in BWR
%    - markers_per_body: the number of points to select in each body for BWR
%    - submatrix_multiplier: the submatrix size for matrix multiplication
%    in BWR is this number times the number of labels. This number must be
%    a divisor of markers_per_body.
%    - d: the dimensionality of the segmentation (2 or 3) (3 default);
%    - debug: print debugging information
% OUTPUT:
%  PR_TABLE: a table containing the various precision-recall curves. The
%  first column is the list of thresholds and the following columns are the
%  precision and recall values for each segmentation at that threshold.
%  COLNAMES: the description of each column in PR_TABLE.

    if ~isfield(spec, 'debug'),
        spec.debug = false;
    end

    % n is the number of segmentations we will be evaluating
    n = numel(spec.ws_filenames);
    if ~isfield(spec, 'ims_mapping'),
        ims_mapping = cell(1, n);
        for i=1:n,
            ims_mapping{i} = 1;
        end
    end
    
    % default values for the index_sets: the full ground truth
    if ~isfield(spec, 'index_sets'),
        spec.index_sets = cell(1, n);
        for i=1:n,
            spec.index_sets{i} = {1:size(spec.gt, 1), 1:size(spec.gt, 2), ...
                1:size(spec.gt, 3)};
        end
    end
    
    % each PR curve may be from a different set of thresholds, but by
    % default, if only one set of thresholds is chosen, it is replicated
    % for each curve
    if ~iscell(spec.ts),
        ts = spec.ts;
        spec.ts = cell(size(spec.ws_filenames));
        for i=1:n,
            spec.ts{i} = ts;
        end
        clear ts;
    end
    ts = spec.ts;
    
    % pr_cell contains 3 columns per precision-recall curve: threshold,
    % precision, recall (in that order). 
    pr_cell = cell(1, 3*n);
    % colnames contains 2 columns per precision-recall curve (precision and
    % recall), plus 1 column at the start for the thresholds (all
    % thresholds are represented in the final pr_table).
    colnames = cell(1, 2*n+1);
    colnames{1} = 'threshold';
    
    % Parameters for body_rand, including defaults
    br_params = struct('perform_watershed_gt', false, ...
                       'perform_watershed_as', false);
    if isfield(spec, 'min_body_size'),
        br_params.min_body_size = spec.min_body_size;
    end
    if isfield(spec, 'markers_per_body'),
        br_params.markers_per_body = spec.markers_per_body;
    end
    if isfield(spec, 'submatrix_multiplier'),
        br_params.submatrix_multiplier = spec.submatrix_multiplier;
    end
    
    % k counts the position in the precision-recall cell array, c counts
    % the position in the colnames array. Although one can be deduced from
    % the other, keeping two counters was clearer and programmatically more
    % convenient
    k=1;
    c=2;
    
    % For each segmentation, and each threshold, compute the body-wise rand
    for i=1:n,
        if spec.debug,
            fprintf('starting evaluation of segmentation %i\n', i);
            t1 = tic;
        end
        m = numel(ts{i});
        pr_cell{k} = ts{i};
        pr_cell{k+1} = zeros([m, 1], 'double');
        pr_cell{k+2} = zeros([m, 1], 'double');
        for j=1:m,
            as = read_segmentation_from_raw(spec.ws_filenames{i}, ...
                sprintf(spec.amb_file_patterns{i}, ts{i}(j)), spec.d);
            [p, r, ~] = body_rand(spec.gt(spec.index_sets{i}{:}), ...
                               as(spec.index_sets{i}{:}), ...
                               br_params);
            pr_cell{k+1}(j) = p;
            pr_cell{k+2}(j) = r;
        end
        colnames{c} = [spec.colname_base{i} '.prec'];
        colnames{c+1} = [spec.colname_base{i} '.rec'];
        k = k+3;
        c = c+2;
        if spec.debug,
            t2 = toc(t1);
            fprintf('done with evaluation %i in %d seconds\n', i, t2);
        end
    end
    pr_table = match_thresholds(pr_cell{:});
    pr_table(:,1) = pr_table(end:(-1):1,1);
end