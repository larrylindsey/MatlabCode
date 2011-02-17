function config = mitochondria_2007KHLateral_s09_1_91(module_id, case_ids)
% config = mitochondria_2007KHLateral_s09_1_91(module_id, case_ids)
% Mitochondria detection on 2007KH Lateral S09 images slices 1 through 91.
%
% module_id
%   1     mitochondria_manual_annotation(config);
%   2     mitochondria_collect_training_samples_intensity_histograms(config);
%   3     Train and apply mitochondria detector using a Leave-one-out
%         protocol. Calls:
%         mitochondria_train_boosted_intensity_hist_detector(config);
%         mitochondria_apply_boosted_intensity_hist_detector(config);
%   4     mitochondria_evaluate_detector(config);
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
% v1  05012008  code borrowed from reconstruction_ex_Alex_volume.m
% v2  05132008  code borrowed from mitochondria_OsO4_1_medulla_ihg.m
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. Constants - Not to be changed by users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config = get_basic_config();
config.DEBUG = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% II. Configuration parameters - to be set by user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Stack parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) stack name
config.stack.name = '2007KHLateral.s09.cropped';
% (b) image name:
% If the files are image_0001.tif, image_0002.tif, etc., then
% the image_prefix would be "image_%04d" and the image_suffix would be
% ".tif".
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
config.stack.image_prefix = 's09%03d.crop';
config.stack.image_suffix = '.tif';
% (c) slice ids: If the images are a001.tif, ..., a010.tif then the
% case_ids are 1,2,...,10.
if(nargin==1)
  config.stack.case_ids = 22:90;
else
  config.stack.case_ids = case_ids;
end
% (d) ROI: Region Of Interest for the stack. If the entire slice image is
% to be processed then assign empty ([]). Otherwise get_roi(config);
config.stack.roi = [];
% config.stack.roi = get_roi(config);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) name for the reconstruction
config.reconstruction.name = 'mitochondria';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Mitochondria
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Whether to print messagse
config.mitochondria.is_verbose = true;
% Whether to display figures
config.mitochondria.is_verbose_figures = true;
% (a) Training image name:
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
config.mitochondria.train.image_prefix = 's09%03d.crop';
config.mitochondria.train.image_suffix = '.tif';
% (b) image ids: : If the images are a006.tif, ..., a010.tif then the
% case_ids are 6,...,10.
config.mitochondria.train.case_ids = 12:21;
% For constructing the feature vector for mitochondria detection.
% () type of features collected for each pixel
config.mitochondria.feature.type = 'heq_intensity_norm_hist_3s2';
% * version name for the feature
config.mitochondria.feature.version = '_v5';
% * The patch size is (2xwindows_size+1) X (2xwindows_size+1).
config.mitochondria.feature.window_sizes = [4 7 10 15 20 25 30];
% * In case of intensity histograms, their bins are specified in intensity_bins
config.mitochondria.feature.intensity_bins = 0:0.1:1;
% * The spatial extent of the Gabor filter (gabor_size+1) X (gabor_size+1).
config.mitochondria.feature.gabor.size = 20;
% * The number theta orientations - response is max for all theta.
config.mitochondria.feature.gabor.n_theta = 8;
% * The spatial frequencies for the Gabor filters
config.mitochondria.feature.gabor.freqs = [0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4];
% () Classifier
% * type
config.mitochondria.model.type = 'boost';
% * version
config.mitochondria.model.version = '_b2';
% * Number of iterations of boosting
config.mitochondria.model.n_iteration = 50;
% * Tree depth
config.mitochondria.model.tree_depth = 3;
% () Suffix for saving mitochondria detection confidences
config.mitochondria.apply.save_suffix = ['.mitochondria_det_conf.hinh3s2', ...
  config.mitochondria.model.version, config.mitochondria.feature.version];
% () Mitochondria detection evaluation
% * version to evaluate
config.mitochondria.evaluate.suffix = config.mitochondria.apply.save_suffix;
% Because of patch-wise features, mitochondria detections are not defined
% with a certain border aroudn the image. This should be the maximum
% window size specified in feature.window_sizes + 5
config.mitochondria.border = 35;
%%%%
% 3.1 Mitochondria segmentation
%%%%
% Method used for segmentation
config.mitochondria.segment.method = 'marker_controlled_watershed';
% Marker Controlled Watershed
% 1. Invert the image - image_i
% 2. Apply imhmin: suppress all local minima in image_i whose depth is less
% than a threshold, minima_depth_threshold (170/255=0.6667): image_i_hmin
% 3. image_i_hmin_i = image_i_hmin - image_i;
% 4. Threshold image_i_hmin_i to prune marker pixels
% Threshold called hmin_threshold. Note that this is relative to
% minima_depth_threshold is not absolute.
%
% Fiter the image before segmentation
config.mitochondria.segment.filter_version = 'heq';
% Suppress minima in the inverted image with depth less than a threshold.
config.mitochondria.segment.minima_depth_threshold = 0.67; % 170/255;
% Apply threshold, hmin_threshold, on inverted image_i_hmin. Pixels with
% values below this threshold are considered to not be markers. Recommended
% value 0.01.
config.mitochondria.segment.hmin_threshold = 0.01;
% Filter sequence for removing very small minima and combinging nearby
% ones. Recommended 'a20_c2_a50'.
config.mitochondria.segment.marker_filter_version = 'a20_c2_a50';
% Mitochondria segment suffix
config.mitochondria.segment.mito_seg_suffix = '.heq.0.67.0.01.a20_c2_a50';
% Fraction area overlap with ground truth threshold. Recommended 0.30
config.mitochondria.segment.gt_area_threshold = 0.30;
% Ground truth suffix
config.mitochondria.segment.gt_suffix = '.0.30';

% 3.2 Classify mitochondria segments
%%%%
% Feature to be used for classification
% E.g., a_cnch_cnch3_va_vna_va3_vna3: area, cumulative normalized
% confidence histogram (1 and 3 sections), vesicle area (pixels) and vesicle
% area normalized (1 and 3 sections).
config.mitochondria.segment.classify.feature.type = 'a_cnch_cnch3_va_vna_va3_vna3';
% Histogram bins for the mitochondria detection values within a segment
config.mitochondria.segment.classify.feature.mito_histogram_bins = -1.0:0.1:1.0;
% Feature version name
config.mitochondria.segment.classify.feature.version = 'v0';
% Mitochondria confidence map to be used in the features.
config.mitochondria.segment.classify.feature.mitochondria_confidence_suffix = ...
  '.hinh3s2_b2_v5';
% Vesicle detection map to be used in the features.
config.mitochondria.segment.classify.feature.vesicle_suffix = '';
% Feature suffix.
config.mitochondria.segment.classify.feature.features_suffix = ...
  '.a_cnch_cnch3_va_vna_va3_vna3.v0.hinh3s2_b2_v5';
% Boosted classifier model parameters. Tree depth, recommended 2
config.mitochondria.segment.classify.model.tree_depth = 2;
% Boosted classifier model parameters. Number of iterations, recommended 2
config.mitochondria.segment.classify.model.n_iteration = 50;
% Boosted classifier model suffix
config.mitochondria.segment.classify.model.model_suffix = '.2.50';
% Segment confidence threshold, recommended 1
config.mitochondria.segment.classify.seg_conf_thresh = -1;
%%%%
% 3.3 Link segments in 3D
%%%%
% Markov Random Field (MRF) parameters. Segment confidence threshold,
% recommended -2.
config.mitochondria.segment.mrf.seg_conf_thresh = -2;
% Markov Random Field parameters. Average linkage weight, recommended 150
config.mitochondria.segment.mrf.avg_linkage_weight = 150;
%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
% Prefixes of suffixes for data generated by different stages of the mitochondria
% detection pipeline.
% E.g., segment labels generated by thresholding the percentage area of
% ground-truth labels would have their suffixes start with
% config.mitochondria.segment.gt_full_suffix. Subsequent stages would
% append config information to this suffix.
%%%%%%%%%%%%%%%%%%%%%%
% For ground truth
config.mitochondria.segment.gt_full_suffix = [...
  '.mito.seg', config.mitochondria.segment.mito_seg_suffix, ...
  '.gt', config.mitochondria.segment.gt_suffix];
% For classifier output
config.mitochondria.segment.classifier_full_suffix = [...
  '.mito.seg', config.mitochondria.segment.mito_seg_suffix, ...
  '.features', config.mitochondria.segment.classify.feature.features_suffix, ...
  '.gt', config.mitochondria.segment.gt_suffix, ...
  '.model', config.mitochondria.segment.classify.model.model_suffix];
% For MRF output
config.mitochondria.segment.mrf_full_suffix = [...
  config.mitochondria.segment.classifier_full_suffix, ...
  '.mrf.', num2str(config.mitochondria.segment.mrf.seg_conf_thresh), ...
  '.', num2str(config.mitochondria.segment.mrf.avg_linkage_weight)];

%%%%%%%%%%%%%%%%%%%%%%
% Suffix given to the label-to-pixelwise-detection routines to generate
% detection maps on the image plane.
%%%%%%%%%%%%%%%%%%%%%%
% If ground-truth labelling is to be employed - useful for baseline.
config.mitochondria.segment.label2detection_suffix = [...
  config.mitochondria.segment.gt_full_suffix];
% For thresholded segment confidences
config.mitochondria.segment.label2detection_suffix = [...
  config.mitochondria.segment.classifier_full_suffix, ...
  '.conf_threshold.', num2str(config.mitochondria.segment.classify.seg_conf_thresh)];
% For MRF output
config.mitochondria.segment.label2detection_suffix = [...
  config.mitochondria.segment.mrf_full_suffix];

%%%%%%%%%%%%%%%%%%%%%%
% Pixelwise detection suffix given to the 3D morphological processing
% routine.
%%%%%%%%%%%%%%%%%%%%%%
config.mitochondria.segment.postprocess_suffix = [...
  config.mitochondria.segment.label2detection_suffix, ...
  '.detections'];

%%%%%%%%%%%%%%%%%%%%%%
% Voxel post processing input suffix
%%%%%%%%%%%%%%%%%%%%%%
% For ground truth, if the 3D morphological processing is _not_ to be
% applied.
config.mitochondria.segment.voxelpostprocess_suffix = [...
  config.mitochondria.segment.gt_full_suffix, ...
  '.detections'];
% For 3D post process output
config.mitochondria.segment.voxelpostprocess_suffix = [...
  config.mitochondria.segment.postprocess_suffix, ...
  '.processed.detections'];

%%%%%%%%%%%%%%%%%%%%%%
% Input suffix for evaluating segmentwise labeling. Semantically this
% stage's input is of the same nature as that of
% config.mitochondria.segment.label2detection_suffix.
%%%%%%%%%%%%%%%%%%%%%%
% For ground truth - should give 100% !
config.mitochondria.segment.evaluate_labeling_suffix = [...
  config.mitochondria.segment.gt_full_suffix];
% For thresholded classifier
config.mitochondria.segment.evaluate_labeling_suffix = [...
  config.mitochondria.segment.classifier_full_suffix, ...
  '.conf_threshold.', num2str(config.mitochondria.segment.classify.seg_conf_thresh)];
% For MRF output
config.mitochondria.segment.evaluate_labeling_suffix = [...
  config.mitochondria.segment.mrf_full_suffix];
% For 3D post process output
config.mitochondria.segment.evaluate_labeling_suffix = [...
  config.mitochondria.segment.mrf_full_suffix, ...
  '.detections.processed'];

%%%%%%%%%%%%%%%%%%%%%%
% Input suffix for evaluating voxelwise detections. Errors within 5
% pixels in the image planes.
%%%%%%%%%%%%%%%%%%%%%%
% For evaluating voxelwise detection before 3D post-processing output
config.mitochondria.segment.evaluate_detection_suffix = [...
  config.mitochondria.segment.label2detection_suffix, ...
  '.detections'];
% For evaluating 3D post-process output
config.mitochondria.segment.evaluate_detection_suffix = [...
  config.mitochondria.segment.voxelpostprocess_suffix, ...
  '.voxelprocessed'];

%%%%%%%%%%%%%%%%%%%%%%
% For displaying detection and 3D volume
%%%%%%%%%%%%%%%%%%%%%%
% For displaying voxelwise detection before 3D post-processing output
config.mitochondria.segment.display_suffix = [...
  config.mitochondria.segment.label2detection_suffix, ...
  '.detections'];
% For displaying 3D post-process output
config.mitochondria.segment.display_suffix = [...
  config.mitochondria.segment.voxelpostprocess_suffix, ...
  '.voxelprocessed'];
% Whether or not to display ground truth also
config.mitochondria.segment.display_gt = false;


%%%%%%%%%%%%%%%%%%%%%%
% For PR curve
%%%%%%%%%%%%%%%%%%%%%%
% For ground truth segment classifications
config.mitochondria.segment.PR_curve_suffixes = [...
  config.mitochondria.segment.evaluate_detection_suffix, ...
  '.PRvoxelresults'];
config.mitochondria.segment.PR_curve_error = '_5';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. Pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pipeline(module_id, config);
beep;

return;
end
