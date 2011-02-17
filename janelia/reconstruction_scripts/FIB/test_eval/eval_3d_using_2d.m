function [precision, recall] = eval_3d_using_2d( ground_truth_mat_fn, ...
    ws_fn, seg_fn, verbose, is_verbose_figures)

    if(nargin<=3)
        is_verbose_figures = false;
    end

    bs = get_3D_stack_from_ws_seg_raw(ws_fn, seg_fn);

    gt = load(ground_truth_mat_fn);

    precision = [0 0];
    recall = [0 0];
    for z=1:length(gt.wcat),
        seg_gt = gt.wcat{z};
        seg_as = bs(:,:,z);
        
        seg_as_b = draw_segment_boundaries_c(double(seg_as));
        seg_as(seg_as_b>0) = 0;
        
        if(is_verbose_figures)
            figure(1);
            greyscale_image = gt.al{z};
            imshow(gt.al{z});
            title('image');
            figure(2);
            imshow(seg_gt, []);
            title('gt');
            figure(3);
            imshow(seg_as, []);
            title('as');
            fprintf('press key to continue\n');
            pause;
        else
            greyscale_image = 0;
        end
        
        % Smooth the segment boundaries to reduce wiggles. Perform watershed one
        % last time to ensure segment boundaries are closed
        b_as = seg_as==0;
        b_as_d = imdilate(b_as, strel('disk', 1));
        seg_as = watershed(b_as_d, 4);
        b_gt = seg_gt==0;
        b_gt_d = imdilate(b_gt, strel('disk', 1));
        seg_gt = watershed(b_gt_d, 4);
        
        % Whether to show match results
        eval_parameters.save_matching = is_verbose_figures;
        % Whether to write match results
        eval_parameters.write_matching = false;
        % Maximum distance (in pixels) upto which pairs of boundary elements are
        % connected in the bipartite graph
        eval_parameters.dmax = 10;
        
        load( ...
            ['/groups/chklovskii/home/nuneziglesiasj/' ...
            'em_reconstructions/reconstructions/spams_denoising/' ...
            'mitochondria/' ...
            sprintf('a.%03d.mitochondria_det_conf.mat', z-1)] ...
        );
        evaluation_mask = filter_image2( ...
            im2double(mitochondria_detection_confidence), 't0.4_a100_d5');
        evaluation_mask = mask_frame(evaluation_mask, 35);
        [pr, rc] = evaluate_segmentation_boundary(...
            seg_as, seg_gt, eval_parameters, im2double(greyscale_image), ...
            evaluation_mask);
        precision = precision + pr;
        recall = recall + rc;
        if verbose,
            fprintf('Slice %d done, current prec: %.2f, current rec: %.2f\n', ...
                z, precision(1)/precision(2), recall(1)/recall(2));
        end

    end

end

