function seg_eval_2007KHApical_s10_subset_gs_amb(module_id)
% seg_eval_2007KHApical_s10_subset_gs_amb(module_id)
% Evaluate segmentation on 2007KHApical.s10.subset_112508 for segmentation
% evaluation paper.
% Algorithm: grayscale agglo. mean boundary
%
% module_id
%   1     mitochondria_manual_annotation(config);
%   2     mitochondria_collect_training_samples_intensity_histograms(config);
%   3     mitochondria_train_boosted_intensity_hist_detector(config);
%   4     mitochondria_apply_boosted_intensity_hist_detector(config);
%   5     segment_2D_grayscale_ladder(config);
%   6     superpixel_2_segment_grayscale_ladder(config);
%   7     linkage_3D_train_intensity_pair_boost(config);
%   8     linkage_3D_gen_linkage_gph_intensity_pair_boost(config);
%   9     generate_al(config);
%         generate_cat(config);
%         generate_superpixel_2_seg_map(config);
%         process_superpixel_2_2Dseg_linkage_graph(config);
%   10    prepare_final_reconstruction_volume(config);
%   11    save_reconstruction_info(config);
%         save_reconstruction_param_to_xml(config);
%   12    get_q_score(config);
%   13    get_rand_score(config);
%   14    get_boundary_match_score(config);
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04382008  Code borrowed from reconstruction_ex_Alex_volume.m
% v1  08052008  code modified for grayscale agglomerative clustering mean
%               boundary
% v2  08252008  apply agglomerative mean boundary segmentation directly on
%               watershed without grayscale-ladder based superpixels
% v3  10312008  seg. eval. for comparitive experiments
% v4  10312008  comparitive experiments - agglo. mean boundary
% v5  11252008  seg. eval. for paper
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
config.stack.name = '2007KHApical.s10.subset_112508';   
% (b) image name:
% If the files are image_0001.tif, image_0002.tif, etc., then
% the image_prefix would be "image_%04d" and the image_suffix would be
% ".tif".
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
config.stack.image_prefix = 'a.%d'; 
config.stack.image_suffix = '.tif';
% (c) slice ids: If the images are a001.tif, ..., a010.tif then the
% case_ids are 1,2,...,10.
config.stack.case_ids = 1:20;
% (d) ROI: Region Of Interest for the stack. If the entire slice image is
% to be processed then assign empty ([]). Otherwise get_roi(config);
config.stack.roi = [];
% config.stack.roi = get_roi(config);
% () Stack resolution parameters for area and boundary histograms. These
% are to be used if the models having been trained on images of resolution
% different from the current one. This might improve the results,
% especially during boot-strapping.
config.stack.length_factor = 1;
config.stack.area_factor = config.stack.length_factor.^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) name for the reconstruction
config.reconstruction.name = 'reconstruction_lamina.OTO.3500x.zhiyuan.May232008.3_5';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Mitochondria
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Training image name:
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
config.mitochondria.train.image_prefix = 'mit_train_redmag_3500_%01d';
config.mitochondria.train.image_suffix = '.tif';
% (b) image ids: : If the images are a006.tif, ..., a010.tif then the
% case_ids are 6,...,10.
config.mitochondria.train.case_ids = 1:1;
% For constructing the feature vector for mitochondria detection.
% () type of features collected for each pixel
config.mitochondria.feature.type = 'heq_intensity_hist_gabor';
% * version name for the feature
config.mitochondria.feature.version = '_v4.1';
% * The patch size is (2xwindows_size+1) X (2xwindows_size+1).
config.mitochondria.feature.window_sizes = [4 10 20 30];
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
config.mitochondria.model.version = '_b';
% * Number of iterations of boosting
config.mitochondria.model.n_iteration = 30;
% * Tree depth
config.mitochondria.model.tree_depth = 3;
% () Suffix for saving mitochondria detection confidences
config.mitochondria.apply.save_suffix = ['.mitochondria_det_conf.hihg', ...
  config.mitochondria.model.version, config.mitochondria.feature.version];
% () Mitochondria detection evaluation
% * version to evaluate
config.mitochondria.evaluate.suffix = config.mitochondria.apply.save_suffix;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Vesicle
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.vesicle.apply.dir = '';
config.vesicle.apply.save_suffix = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Superpixel - multi stage
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel.method = 'grayscale_ladder';
% (b) thresholds on the boundary field
config.superpixel.f_thresholds = 0.2:0.05:0.9;
% (c) minimum area of a segment - below this they merged with neighbors
config.superpixel.area_thresholds = 400;
% (d) Whether mitochondria are to be used to suppress false boundaries
config.superpixel.use_mitochondria = false;
% (e) The confidence threshold to be applied on the mitochondria
config.superpixel.mitochondria_confidence_threshold = 0.0; 
% (f) The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel.mitochondria_erosion = 5; % 20
% * Minimum of a thresholded connected component for it to be considered
% mitochondria
config.superpixel.mitochondria_min_area = 40;
% (g) Whether to use vesicle detection to suppress false boundaries
config.superpixel.use_vesicle = false;
% (h) The threshold for the vesicle detection confidence
config.superpixel.vesicle_threshold = 0;
% (i) Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.
config.superpixel.vesicle_dilation = 0;
% (j) Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel.filter_version = 'v_heq_mf2';
% () whether to display intermediate results
config.superpixel.is_verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Superpixel to segment
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel_2_seg.method = 'grayscale_AglMeanB';
% () Superpixel method
config.superpixel_2_seg.superpixel_method = 'grayscale_ladder';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria m<threshold_d<amount of erosion>>
config.superpixel_2_seg.superpixel_suffix = '.gs_l_T0.55_L200.v_heq_mf2'; 
% (c) thresholds on the boundary field
config.superpixel_2_seg.f_threshold_seq = 0.55:0.01:0.72;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel_2_seg.save_f_thresholds = 0.62:0.02:0.72;
% () minimum boundary length (pixel units) for two segments to be merged
% should increase with increasing f_thresholds. E.g., 0 for f_thresholds
% <=0.65 and 20 otherwise
config.superpixel_2_seg.length_threshold_seq = ...
  25*(config.superpixel_2_seg.f_threshold_seq>0.65);
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel_2_seg.use_mitochondria = false;
% () The confidence threshold to be applied on the mitochondria
config.superpixel_2_seg.mitochondria_confidence_threshold = 0.0;
% () The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel_2_seg.mitochondria_erosion = 20;
% () Whether to use vesicle detection to suppress false boundaries
config.superpixel_2_seg.use_vesicle = false;
% () The threshold for the vesicle detection confidence
config.superpixel_2_seg.vesicle_threshold = 0;
% () Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.
config.superpixel_2_seg.vesicle_dilation = 0;
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel_2_seg.filter_version = 'v_heq_mf2';
% () whether to display intermediate results
config.superpixel_2_seg.is_verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Evaluate segment boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%
% () Ground truth segmentation to be used for evaluation.
% Format similar to that for volume reconstruction evaluation.
% seg is a 3D matrix of body labels for the voxels.
config.eval_segment_boundary.groundtruth_file = 'proofread.manual.seg.mat';
% () Lsit of segments in the ground-truth that may be over-segmented without
% penalty. Use this for irrelevant structurs , e.g., glia.
% config.eval_segment_boundary.oversegment_ignore_groundtruth_label = [];
% () Method to be evaluated
config.eval_segment_boundary.method = 'grayscale_AglMeanB';
% () Parameters for the method to be evaluated in the segmentation suffix
% There are two options:
% * Verbose suffix to provide a specific combination of parameters
% config.eval_segment_boundary.seg_suffix = '.gs_l_T0.56_L600.v_heq_mf10_m0';
% * Or, specify a list of parameters. In this case, provide the ANSI C
% compatible printf format string and the list of ids to be printed with
% the format string.
config.eval_segment_boundary.seg_suffix_format = ...
  ['.gs_amb_sp_T%s', config.superpixel_2_seg.superpixel_suffix, '.v_heq_mf2'];
s = {};
for f = config.superpixel_2_seg.save_f_thresholds
  s{end+1} = [num2str(f,'%g'),'_',num2str(f-0.01,'%g'),'_b',num2str(25*(f>0.65),'%d')];
end
config.eval_segment_boundary.seg_suffix_id =  s;
% () Maximum matching distance between ground-truth and automatic
% segmentation boundaries. Depends upon resolution. Keeping it very small
% results in zig-zag boundaries being penalized even though the overall
% segmentation is OK.
config.eval_segment_boundary.max_match_dist = 15;
% () Verbose
config.eval_segment_boundary.is_verbose = true;
% () Mode of evaluation: 'pr-curve', 'match-save'
% pr-curve: A precision-recall is generated, typically for a range of
% segmenation parameters specified using seg_suffix_id.
% match-save: Saves the bipartite matching results to .tif images and
% creates a .tex file. Compiling this file using pdflatex produces a pdf of
% the matching results. This enables detailed inspection of the evaluation.
config.eval_segment_boundary.mode = 'pr-curve';
% config.eval_segment_boundary.mode = 'match-save';
% type of marker for the pr-curve
config.eval_segment_boundary.marker = 'or-';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. Reconstruction Routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pipeline(module_id, config);

return

end