image = imread('/groups/chklovskii/medulla/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/a.012.tif');
interior_map = imread('/groups/chklovskii/chklovskiilab/RandomForests/classification result2.tif');
mito = load('/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/mitochondria.shiv/a.012.mitochondria_det_conf.mat');

[label_map, watershed_label_map] = get_segment_map_from_random_forest_interior_map(...
  image, interior_map, ...
  'mitochondria_map', mito.mitochondria_detection_confidence>0, ...
  'interior_area_threshold', 0);

seg_gt = load('/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/2D_segmentation_results/proofread_result/a.012.prfrd_seg.satoko.031009.mat');
plot_segment_boundaries(image, seg_gt.label_map, 1, 8);
title('ground-truth segment map boundaries');








fprintf('Evaluating segment map using bipartite matching\n');
% maximum distance between boundary points for them to be allowed to match
% [pixel units]
evaluate_config.dmax = 15; 
evaluate_config.is_verbose = true;
evaluate_config.is_verbose_figure = false;
boundary_match_result = evaluate_boundary_bipartite_match(...
  label_map==0, seg_gt.label_map==0, evaluate_config, zeros(size(label_map)));

%%%%%%%%%%%%%% Results %%%%%%%%%%%%%%
%%% a.012.tif
%%%   classification result.tif   recall = 0.9769, precision = 0.8856    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








fprintf('Evaluating segment boundaries using pixel overlap (dilated)\n');
dilation_radius = 15;
gt_boundaries = double(seg_gt.label_map==0);
auto_boundaries = double(watershed_label_map==0);
gt_boundaries_detected = gt_boundaries .* ...
  imdilate(auto_boundaries, strel('disk', dilation_radius));
auto_boundaries_detected = auto_boundaries .* ...
  imdilate(gt_boundaries, strel('disk', dilation_radius));
gt_boundaries_missed = gt_boundaries .* (1 - gt_boundaries_detected);
gt_boundaries_missed = bwareaopen(gt_boundaries_missed, 2);
plot_segment_boundaries(image, 1 - gt_boundaries_missed, 1, 9);
title('missed ground-truth boundaries in watershed');
auto_boundaries_false = auto_boundaries .* (1 - auto_boundaries_detected);
auto_boundaries_false = bwareaopen(auto_boundaries_false, 2);
plot_segment_boundaries(image, 1 - auto_boundaries_false, 1, 10);
title('false watershed segment boundaries');

recall = nnz(gt_boundaries_detected) / nnz(gt_boundaries);
precision = nnz(auto_boundaries_detected) / nnz(auto_boundaries);
fprintf('watershed: recall: %g, precision: %g\n', recall, precision);

fprintf('Evaluating segment boundaries using pixel overlap (dilated)\n');
dilation_radius = 15;
gt_boundaries = double(seg_gt.label_map==0);
auto_boundaries = double(label_map==0);
gt_boundaries_detected = gt_boundaries .* ...
  imdilate(auto_boundaries, strel('disk', dilation_radius));
auto_boundaries_detected = auto_boundaries .* ...
  imdilate(gt_boundaries, strel('disk', dilation_radius));
gt_boundaries_missed = gt_boundaries .* (1 - gt_boundaries_detected);
gt_boundaries_missed = bwareaopen(gt_boundaries_missed, 2);
plot_segment_boundaries(image, 1 - gt_boundaries_missed, 1, 11);
title('missed ground-truth boundaries in ladder');
auto_boundaries_false = auto_boundaries .* (1 - auto_boundaries_detected);
auto_boundaries_false = bwareaopen(auto_boundaries_false, 2);
plot_segment_boundaries(image, 1 - auto_boundaries_false, 1, 12);
title('false ladder segment boundaries');

recall = nnz(gt_boundaries_detected) / nnz(gt_boundaries);
precision = nnz(auto_boundaries_detected) / nnz(auto_boundaries);
fprintf('ladder: recall: %g, precision: %g\n', recall, precision);

%%%%%%%%%%%%%% Results %%%%%%%%%%%%%%
%%% a.012.tif
%%%   classification result.tif   recall = 0.9815, precision = 0.9398
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
