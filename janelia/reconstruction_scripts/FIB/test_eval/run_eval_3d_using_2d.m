function run_eval_3d_using_2d( params, varargin )
% Function RUN_EVAL_3D_USING_2D( PARAMS, VARARGIN )
% Run the evaluation pipeline on a number of stacks (varargin are the stack
% IDs.)
%  - INPUT:
%      * params: the parameters for the evaluation. These include:
%          - gt_mat_fn: the ground truth wcat.mat file (e.g.
%          proofread.victor.wcat.mat)
%          - base_dir: the directory in which the segmentation results are
%          found
%          - watershed_base: where inside the base directory the watershed
%          results are found (subdirectory + filename prefix)
%          - segment_base: as above, but for superpixel-segment map
%          - f_threshold_ambcs_list: the list of thresholds used for ambcs
%          - dim: The dimension of the signal denoising
%          - verbose: print out runtime messages
%          - verbose_figs: print out figures of the result
%      * varargin: the stack ids to be processed
%  - OUTPUT: None.
%  - SIDE EFFECTS: Write out evaluation statistics to disk

    if isdeployed,
        params = eval(params);
    end

    if isnumeric(params),
        params = struct();
    end

    if ~isfield(params, 'gt_mat_fn'),
        gt_mat_fn = '/groups/chklovskii/home/vitaladevunis/temp/saundersm.tiledtestfib2/proofread.saundersm.wcat.mat';
    else
        gt_mat_fn = params.gt_mat_fn;
    end

    if ~isfield(params, 'base_dir'),
        base_dir = '/groups/chklovskii/home/nuneziglesiasj/em_reconstructions/reconstructions/spams_denoising/3D_segmentation_results';
    else
        base_dir = params.base_dir;
    end

    if ~isfield(params, 'watershed_base'),
        watershed_base = 'watershed/seg_stack.0.203';
    else
        watershed_base = params.watershed_base;
    end

    if ~isfield(params, 'segment_base'),
        segment_base = 'agglo_mean_boundary/seg_stack.0.203';
    else
        segment_base = params.segment_base;
    end

    if ~isfield(params, 'f_threshold_ambcs_list'),
        f_threshold_ambcs_list = 70:210;
    else
        f_threshold_ambcs_list = params.f_threshold_ambcs_list;
    end

    if ~isfield(params, 'dim'),
        params.dim = 3;
    end

    if ~isfield(params, 'verbose'),
        params.verbose = 0;
    end
    
    if ~isfield(params, 'verbose_figs'),
        params.verbose_figs = 0;
    end

    stack_id_list = varargin;

    result_matrix = zeros( ...
                        numel(stack_id_list),numel(f_threshold_ambcs_list), 4);

    t0 = tic;
    for j=1:numel(stack_id_list),
        t1 = tic;
        ws_fn = [base_dir '/' watershed_base '.' stack_id_list{j} '.ws.T10.raw'];
        for k=1:numel(f_threshold_ambcs_list),
            f_threshold_ambcs = f_threshold_ambcs_list(k);
            seg_fn = sprintf([base_dir '/' segment_base '.' stack_id_list{j} ...
                '.ws.T10.ambl.T90.L50.A1000.ambc.T%d.0_203.raw'], ...
                f_threshold_ambcs);
            if params.verbose,
                fprintf('filenames:\n   %s\n   %s\n   %s\n', ...
                                                gt_mat_fn, ws_fn, seg_fn);
            end
            [precision, recall] = eval_3d_using_2d(gt_mat_fn, ws_fn, ...
                              seg_fn, params.verbose, params.verbose_figs);
            result_matrix(j,k,:) = [precision, recall];
        end
        t2 = toc(t1);
        tt = toc(t0);
        for k=1:size(result_matrix,3),
        dlmwrite(['~/Projects/em_denoising/data/precision-recall-2D.' ...
            stack_id_list{1} '.' ...
            num2str(f_threshold_ambcs_list(1)) '.' num2str(k) '.txt'], ...
            result_matrix(:,:,k));
        end
        if params.verbose,
            fprintf(['done with stack %s. Time taken: %.2f.' ...
                ' Total time: %.2f.\n'], stack_id_list{j}, t2/3600, tt/3600);
        end
    end

end
