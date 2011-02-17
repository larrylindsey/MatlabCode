function [pr_table, colnames] = all_cepr(spec)
% [PR_TABLE, COLNAMES] = ALL_CEPR(SPEC)
% Generate a table of close-enough precision-recall values encompassing 
% multiple segmentation methods, potentially from different parts of an
% image stack.
% INPUT:
%  SPEC: a structure with the following fields. 
%    - ws_filenames: cell array of the filenames for the watershed files
%    - amb_file_patterns: cell array of the amb files corresponding to the
%    watershed files
%    - ts: cell array of the thresholds for each amb_file_pattern OR single
%    set of thresholds, identical for each amb_file_pattern
%    - gt: ground truth boundary map
%    - se: the structuring elements to do the evaluation dilation
%    - se_mapping: which structuring elements should be used for each
%    ws_filename
%    - ims: cell array of images for the visualization
%    - ims_mapping: which image stack from ims to use for each
%    visualization
%    - index_sets: which parts of the image to evaluate
%    - vizfns: the base filename for each visualization stack
%    - vizdir: the directory in which to put the visualization stacks
%    - colname_base: the base name of each column
%    - trim: the amount of margin to remove from evaluation in CEPR
% OUTPUT:
%  PR_TABLE: a table containing the various precision-recall curves. The
%  first column is the list of thresholds and the following columns are the
%  precision and recall values for each segmentation at that threshold.
%  COLNAMES: the description of each column in PR_TABLE.

    % n is the number of segmentations we will be evaluating
    n = numel(spec.ws_filenames);
    if ~isfield(spec, 'ims_mapping'),
        ims_mapping = cell(1, n);
        for i=1:n,
            ims_mapping{i} = 1;
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
    
    % m is the number of PR curves we will be "drawing". This is the sum
    % over all segmentations of the number of different structuring
    % elements used to evaluate them
    m = 0;
    for i=1:numel(spec.se_mapping),
        m = m+numel(spec.se_mapping{i});
    end
    
    % pr_cell contains 3 columns per precision-recall curve: threshold,
    % precision, recall (in that order). 
    pr_cell = cell(1, 3*m);
    % colnames contains 2 columns per precision-recall curve (precision and
    % recall), plus 1 column at the start for the thresholds (all
    % thresholds are represented in the final pr_table).
    colnames = cell(1, 2*m+1);
    colnames{1} = 'threshold';
    
    % k counts the position in the precision-recall cell array, c counts
    % the position in the colnames array. Although one can be deduced from
    % the other, keeping two counters was clearer and programmatically more
    % convenient
    k=1;
    c=2;
    
    % For each segmentation, compute the boundary confidence map, then for
    % each structuring element, compute the close-enough precision recall
    % (CEPR) curve, and finally join all the precision-recall curves under
    % a common scale
    for i=1:n,
        as = tsegments2beckham(spec.ws_filenames{i}, ...
                               spec.amb_file_patterns{i}, ...
                               spec.ts{i} ...
                               );
        for j=1:numel(spec.se_mapping{i}),
            cur_se = spec.se{spec.se_mapping{i}(j)};
            if numel(spec.ims{spec.ims_mapping{i}}) > 1,
                cur_im = spec.ims{spec.ims_mapping{i}}(spec.index_sets{i}{:});
            else
                cur_im = 0;
            end
            [p, r, t, v] = close_enough_precision_recall( ...
                spec.gt(spec.index_sets{i}{:}), ...
                as(spec.index_sets{i}{:}), 1, cur_se, cur_im, spec.trim);
            pr_cell{k} = t;
            pr_cell{k+1} = p;
            pr_cell{k+2} = r;
            k = k+3;
            colnames{c} = sprintf([spec.colname_base{i} '.prec.se%i'], ...
                spec.se_mapping{i}(j));
            colnames{c+1} = sprintf([spec.colname_base{i} '.rec.se%i'], ...
                spec.se_mapping{i}(j));
            c = c+2;
            if numel(v) > 1,
                for l=1:size(v, 5),
                    eval(['!mkdir -p ' spec.vizdir '/' spec.vizfns{i}]);
                    write_image_stack(squeeze(v(:,:,:,:,l)), struct( ...
                        'filename', ...
                        sprintf([spec.vizfns{i} '.se%i.t%03i.tif'], j, t(l)), ...
                        'directory', [spec.vizdir '/' spec.vizfns{i}]));
                end
            end
        end
    end
    pr_table = match_thresholds(pr_cell{:});
end
