function [prec, rec, ts, viz] = close_enough_precision_recall(gt, as, ...
                                    high_boundary, r, viz, trim)
% [PREC, REC, TS, VIZ] = CLOSE_ENOUGH_PRECISION_RECALL(GT, PRED [,  
%                                               HIGH_BOUNDARY, R, VIZ, TRIM])
% Compute precision and recall by approximately matching boundaries using
% dilation.
% INPUT:
%   - GT: The ground truth volume boundary map
%   - PRED: The predicted boundary map. This map can contain continuous
%   predictive values, such as mean-at-merge values, for the boundary voxels.
%   - HIGH_BOUNDARY: defines both 1) the label of boundary pixels in GT and,
%   2) that high values in PRED are more likely to be boundary than low values.
%   - R: either a radius for a spherical 3D structuring
%   element, or an actual structuring element, for the dilation operation.
%   - VIZ: either 0 (for no visualization) or a grayscale volume (output a
%   4D visualization of the volume) (TP, FP and FN boundaries for each
%   threshold used in the precision-recall).
%   - TRIM: amount of margin to remove from each face of the volume when
%   performing evaluation.

    if nargin < 3,
        high_boundary = 1;
    end
    if nargin < 4,
        r = zeros([3 3 3], 'uint8');
        r(:,:,2) = 1;
        r(:,2,:) = 1;
        r(2,:,:) = 1;
    end
    if nargin < 5,
        viz = 0;
    end
    if nargin < 6,
        trim = 10;
    end
    
    s = size(gt);

    if trim > 0,
        gt = gt(1+trim:s(1)-trim, 1+trim:s(2)-trim, 1+trim:s(3)-trim);
        as = as(1+trim:s(1)-trim, 1+trim:s(2)-trim, 1+trim:s(3)-trim);
        if numel(viz) > 1,
            viz = viz(1+trim:s(1)-trim, 1+trim:s(2)-trim, 1+trim:s(3)-trim);
        end
        s = size(gt);
    end
    
    % change to high predictive values mean boundary if it is not the case
    if ~high_boundary,
        as = max(as(:))-as;
    end
    
    % if 'r' is a scalar, make a ball structuring element.
    % Otherwise, use 'r' directly as a structuring element.
    if numel(r) == 1 && ~strcmp(class(r), 'strel'),
        se = strel('ball', r, 2*r+1);
    else
        se = r;
    end
    
    % the method depends on expanding the region matching a boundary. If
    % the boundary label is 0 in gt, then we should do an erosion
    % instead of a dilation. However, it should always be a dilation in the
    % predictor volume, since we reversed the sign of the predictor in case
    % of high_boundary == 0;
    if ~high_boundary,
        gtd = imerode(gt, se);
    else
        gtd = imdilate(gt, se);
    end
    asd = imdilate(as, se);
    
    % IMPORTANT: the method requires that 
    %            unique(as) == unique(asd)
    % this statement ensures that, while preserving the conceptual
    % correctness of the algorithm.
    asd(logical(as)) = as(logical(as));
    ts = sort(unique(as)); ts = ts(2:end);
    
    [~, prec] = perfcurve(gtd(:), double(as(:)), ...
                    high_boundary, 'xCrit', 'reca', 'yCrit', 'prec');
    [rec, ~] = perfcurve(gt(:), double(asd(:)), ...
                    high_boundary, 'xCrit', 'reca', 'yCrit', 'prec');
                
    % perfcurve returns the trivial (recall, precision) points (0,1) and
    % (1,0), except (0,1) is (0,NaN). These are annoying and so are tossed
    % out by this function.
    prec = prec(2:end-1);
    rec = rec(2:end-1);
    
    if numel(viz) > 1,
        viz = reshape(viz, [s(1), s(2), 1, s(3) 1]);
        viz_orig = viz;
        viz = repmat(viz, [1 1 3 1 numel(ts)]);
        h = waitbar(0, 'Building 4D precision-recall visualization...');
        for i=1:numel(ts),
            viztp = viz_orig;
            vizfp = viz_orig;
            vizfn = viz_orig;
            as_cur = uint8(as >= ts(i));
            asd_cur = uint8(asd >= ts(i));
            tp = reshape(as_cur .* gtd + gt .* asd_cur, [s(1) s(2) 1 s(3) 1]);
            tp(logical(tp)) = 1;
            fp = reshape(as_cur .* (1-gtd), [s(1) s(2) 1 s(3) 1]);
            fn = reshape(gt .* (1-asd_cur), [s(1) s(2) 1 s(3) 1]);
            fp = fp .* (1-tp);
            fn = fn .* (1-tp);
            vox = logical(tp+fp+fn);
            viztp(vox) = 255*tp(vox);
            vizfp(vox) = 255*fp(vox);
            vizfn(vox) = 255*fn(vox);
            viz(:,:,1,:,i) = vizfn; % fn in red
            viz(:,:,2,:,i) = vizfp; % fp in green
            viz(:,:,3,:,i) = viztp; % tp in blue
            waitbar(i/numel(ts), h);
        end
        delete(h);
    end
end