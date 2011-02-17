function seg_eval_lamina_leginon_cropA_03262009_gs_amb_l2(module_id)
% seg_eval_lamina_leginon_cropA_03262009_gs_amb_l(module_id)
% Evaluate segmentation on lamina.leginon.cropA.03262009 stack
% Algorithm: grayscale ladder
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
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. Constants - Not to be changed by users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config = get_basic_config();
global config_global
% Whether to use hashing when generating file names during storage and
% retrieval. This reduces the filename lengths.
config_global.storage.file_name_is_hashed__ = false;

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
config.stack.name = 'lamina.leginon.cropA.03262009';   
% (b) image name:
% If the files are image_0001.tif, image_0002.tif, etc., then
% the image_prefix would be "image_%04d" and the image_suffix would be
% ".tif".
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
config.stack.image_prefix = 'a.%03d'; 
config.stack.image_suffix = '.tif';
% (c) slice ids: If the images are a001.tif, ..., a010.tif then the
% case_ids are 1,2,...,10.
config.stack.case_ids = 1:3;
% (d) ROI: Region Of Interest for the stack. If the entire slice image is
% to be processed then assign empty ([]). Otherwise get_roi(config);
config.stack.roi = [];
% config.stack.roi = get_roi(config);
% () Stack resolution parameters for area and boundary histograms. These
% are to be used if the models having been trained on images of resolution
% different from the current one. This might improve the results,
% especially during boot-strapping.
% config.stack.length_factor = 1;
% config.stack.area_factor = config.stack.length_factor.^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) name for the reconstruction
config.reconstruction.name = 'reconstruction_lamina.OTO.3500x.zhiyuan.May232008.3_5';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Mitochondria
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Training image name:
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
config.mitochondria.train.image_prefix = 'a%03d';
config.mitochondria.train.image_suffix = '.tif';
% (b) image ids: : If the images are a006.tif, ..., a010.tif then the
% case_ids are 6,...,10.
config.mitochondria.train.case_ids = 1:3;
% (c) For constructing the feature vector for mitochondria detection. In
% case intensity histograms, the patch size is (2xwindows_size+1) X (2xwindows_size+1).
% The intensity histogram's bins are specified in intensity_bins
config.mitochondria.feature.type = 'heq_intensity_hist';
config.mitochondria.feature.window_sizes = [15 25 35];
config.mitochondria.feature.intensity_bins = 0:0.1:1;
% (d) Classifier type
config.mitochondria.model.type = 'boost';
% (e) Number of iterations of boosting
config.mitochondria.model.n_iteration = 30;
% (f) Tree depth
config.mitochondria.model.tree_depth = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Vesicle
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.vesicle.apply.dir = '';
config.vesicle.apply.save_suffix = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Superpixel
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.superpixel = repmat(config.superpixel, [2 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.1 Superpixel - agglomerative mean boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel(1).method = 'grayscale_AglMeanB';
% () thresholds on the boundary field
config.superpixel(1).f_threshold_seq = 0.1:0.01:0.66; %0.9;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel(1).save_f_thresholds = 0.4:0.02:0.6;
% () minimum boundary length (pixel units) for two segments to be merged
% should increase with increasing f_thresholds. E.g., 0 for f_thresholds
% <=0.65 and 20 otherwise
config.superpixel(1).length_threshold_seq = ...
  25*(config.superpixel(1).f_threshold_seq>0.65);
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel(1).use_mitochondria = true;
% () The confidence threshold to be applied on the mitochondria
config.superpixel(1).mitochondria_confidence_threshold = 0.5; 
% () The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel(1).mitochondria_erosion =20; % 20
% * Minimum of a thresholded connected component for it to be considered
% mitochondria
config.superpixel(1).mitochondria_min_area = 40;
% () Whether to use vesicle detection to suppress false boundaries
config.superpixel(1).use_vesicle = false;
% () The threshold for the vesicle detection confidence
config.superpixel(1).vesicle_threshold = 0;
% () Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.
config.superpixel(1).vesicle_dilation = 0;
% % Image filtering for boundary detection to be used for watershed
% % boundary = 1 - filtered_image
% config.superpixel(1).watershed_filter_version = 'heq2_o1';
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel(1).filter_version = 'v_heq2_mf2_e2';
% () whether to display intermediate results
config.superpixel(1).is_verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.2 Superpixel - grayscale ladder
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel(2).method = 'grayscale_ladder';
% () Superpixel method
config.superpixel(2).superpixel_method = 'grayscale_AglMeanB';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
s = {};
for f = config.superpixel(1).save_f_thresholds
  s{end+1} = ['.gs_amb_T', num2str(f,'%g'), '_', num2str(f-0.01,'%g'), ...
    '_b', num2str(25*(f>0.65),'%d'), '.v_heq2_mf2_e2_m0.5_d20_a40.v_heq2_mf2_e2'];
end
config.superpixel(2).superpixel_suffix = s; 
% () thresholds on the boundary field
config.superpixel(2).f_thresholds = 0.2:0.01:0.5;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel(2).save_f_thresholds = config.superpixel(2).f_thresholds;
% () minimum area of a segment - below this they merged with neighbors
config.superpixel(2).area_thresholds = 150;
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel(2).use_mitochondria = true;
% () The confidence threshold to be applied on the mitochondria
config.superpixel(2).mitochondria_confidence_threshold = 0.5; 
% () The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel(2).mitochondria_erosion = 20; % 20
% * Minimum of a thresholded connected component for it to be considered
% mitochondria
config.superpixel(2).mitochondria_min_area = 40;
% () Whether to use vesicle detection to suppress false boundaries
config.superpixel(2).use_vesicle = false;
% () The threshold for the vesicle detection confidence
config.superpixel(2).vesicle_threshold = 0;
% () Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.
config.superpixel(2).vesicle_dilation = 0;
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel(2).filter_version = 'v_heq2_mf2_e2';
% () whether to display intermediate results
config.superpixel(2).is_verbose = true;

for i=0.4:0.02:0.6
  i
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Evaluate segment boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%
% () Ground truth segmentation to be used for evaluation.
% Format similar to that for volume reconstruction evaluation.
% seg is a 3D matrix of body labels for the voxels.
config.eval_segment_boundary.groundtruth_file = 'marta.proofread.20090213.seg.mat';
% () Lsit of segments in the ground-truth that may be over-segmented without
% penalty. Use this for irrelevant structurs , e.g., glia.
config.eval_segment_boundary.oversegment_ignore_groundtruth_label = 74;
% () Method to be evaluated
% () Parameters for the method to be evaluated in the segmentation suffix
% There are two options:
% * Verbose suffix to provide a specific combination of parameters
% config.eval_segment_boundary.seg_suffix = '.gs_l_T0.56_L600.v_heq_mf10_m0';
% * Or, specify a list of parameters. In this case, provide the ANSI C
% compatible printf format string and the list of ids to be printed with
% the format string.
% config.eval_segment_boundary.method = 'grayscale_AglMeanB';
% config.eval_segment_boundary.seg_suffix_format = ...
%   [ '.gs_amb%s', ...
%     '.v_heq2_mf2', ...
%     '.v_heq2_mf2'];
config.eval_segment_boundary.method = 'grayscale_ladder';
config.eval_segment_boundary.seg_suffix_format = ...
  [ '.gs_l_sp%s', ...
    '.v_heq2_mf2_e2_m0.5_d20_a40', ...
    '.v_heq2_mf2_e2', ...
    '.v_heq2_mf2_e2_m0.5_d20_a40'];
s = {};
f1 = i
for f2 = config.superpixel(2).save_f_thresholds
  s{end+1} = ['_T', num2str(f2, '%g'), '_L150', ...
    '.gs_amb_T', num2str(f1,'%g'), '_', num2str(f1-0.01,'%g'), ...
    '_b', num2str(25*(f1>0.65),'%d')];
end
config.eval_segment_boundary.seg_suffix_id = s;
% () Maximum matching distance between ground-truth and automatic
% segmentation boundaries. Depends upon resolution. Keeping it very small
% results in zig-zag boundaries being penalized even though the overall
% segmentation is OK.
config.eval_segment_boundary.max_match_dist = 15;
% () Verbose
config.eval_segment_boundary.is_verbose = false;
% () Mode of evaluation: 'pr-curve', 'match-save'
% pr-curve: A precision-recall is generated, typically for a range of
% segmenation parameters specified using seg_suffix_id.
% match-save: Saves the bipartite matching results to .tif images and
% creates a .tex file. Compiling this file using pdflatex produces a pdf of
% the matching results. This enables detailed inspection of the evaluation.
 config.eval_segment_boundary.mode = 'pr-curve';
%config.eval_segment_boundary.mode = 'match-save';
% type of marker for the pr-curve
config.eval_segment_boundary.marker = 'or-';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. Reconstruction Routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pipeline(module_id, config);


end
return
end