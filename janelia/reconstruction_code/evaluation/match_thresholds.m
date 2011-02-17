function pr_table = match_thresholds(varargin)
% PR_TABLE = MATCH_THRESHOLDS(VARARGIN)
% Create a unified table containing a single list of thresholds from any
% number of (thresholds, precision, recall) column triplets, given as a
% variable-length argument list, in that order.

    vs = varargin; % for notational convenience only
    
    % compute the union of all sets of thresholds for the unified table
    all_thresholds = [];
    for i=1:(numel(vs)/3),
        all_thresholds = union(all_thresholds, vs{3*i-2});
    end
    
    % compute precision recall at *all* thresholds. PR at an un-represented
    % threshold is equal to PR at the next-highest represented threshold.
    % If the highest threshold is not represented, we set precision to be 0
    % and recall to be 1. 
    pr_table = zeros([numel(all_thresholds) (numel(varargin)*2/3+1)]);
    pr_table(:,1) = all_thresholds;
    for i=1:(numel(vs)/3),
        p = (-1)*ones([numel(all_thresholds) 1]);
        r = p;
        for j=1:numel(p),
            found = find(vs{3*i-2} == all_thresholds(j), 1, 'first');
            if ~isempty(found),
                p(j) = vs{3*i-1}(found);
                r(j) = vs{3*i}(found);
            end
        end
        cur_p = 1;
        cur_r = 0;
        for j=numel(p):(-1):1,
            if p(j) == -1,
                p(j) = cur_p;
                r(j) = cur_r;
            else
                cur_p = p(j);
                cur_r = r(j);
            end
        end
        vs{3*i-1} = p;
        vs{3*i} = r;
    end
    for i=1:(numel(vs)/3),
        pr_table(:,2*i) = vs{3*i-1};
        pr_table(:,2*i+1) = vs{3*i};
    end
    
    % The lowest recall is actually at the highest threshold, so we invert
    % the threshold list in the final table.
    pr_table(:,1) = pr_table(end:-1:1, 1);
end
