function config = seg_eval_medulla_Leginon_misc_LDA_v0(config, case_ids, is_verbose, is_verbose_figures)
% seg_eval_medulla_Leginon_3500x_fall08_misc_BEL_amb_l_abifl(module_id)
% Evaluate LDA for segmentation boundary on misc. dataset from 
% medulla.HPF.Leginon.3500x.zhiyuan.fall2008.
% Algo: LDA + agglo. mean boundary + ladder + ABIFL
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
% v0  04162008  init code
% v1  05062008  code borrowed from reconstruction script
%               Marta's reconstruction_lamina_OTO_Jan222008_ministack.m
% v2  05072008  Modified for segmentation evaluation
% v3  06042008  Modified for medulla HPF FS
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. Constants - Not to be changed by users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Whether to use hashing when generating file names during storage and
% retrieval. This reduces the filename lengths.
global config_global
config_global.storage.file_name_is_hashed__ = false;

config.stack.dir = '/groups/chklovskii/medulla/';

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
config.stack.name = 'medulla.HPF.Leginon.3500x.zhiyuan.fall2008';   
% (b) image name:
% If the files are image_0001.tif, image_0002.tif, etc., then
% the image_prefix would be "image_%04d" and the image_suffix would be
% ".tif".
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
% config.stack.image_prefix = 'a.012.cropped%03d';  % 'a.%03d';
config.stack.image_prefix = 'a.%03d';  % 'a.%03d';  

config.stack.image_suffix = '.tif';
% (c) slice ids: If the images are a001.tif, ..., a010.tif then the
% case_ids are 1,2,...,10.
if(nargin<2 || isempty(case_ids))
  config.stack.case_ids = 12;
else
  config.stack.case_ids = case_ids;
end
% Option 2: If the images correspond to a trakEM project then
% * specify the trakEM file structure xml file name
% % % config.stack.image_structure = 'crop4_global_alignment_0161_0860.xml';
% % % % * slice ids: specify the z-plane numbers.
% % % config.stack.case_ids = 443;
% (d) ROI: Region Of Interest for the stack. If the entire slice image is
% to be processed then assign empty ([]). Otherwise get_roi(config);
config.stack.roi = [];
% config.stack.roi.xmin = 101;
% config.stack.roi.xmax = 1800;
% config.stack.roi.ymin = 101;
% config.stack.roi.ymax = 1800;
% config.stack.roi = get_roi(config);
% (e) Alignment parameters to describe the registration of different slices
% with respect to each other. The alignment should be done in IMOD or
% equivalent program to generate a .xf file.
% A stack is prealigned if the images in the stack's directory are
% already aligned to each other. If the stack is not prealigned then the
% 2D segmentation (superpixel and superpixel-to-segments) is performed on
% the unaligned images. The alignment parameters are used for linkage and
% to construct the proofreading data structures.
% * If the stack is prealigned then set is_prealigned as true else false
config.stack.align.is_prealigned = true;
% * If not prealigned then name the .xf file specifying the alignment.
% Assumed to be located in the stack directory.
% config.stack.align.xf_file_name = '1.xf';
% config.stack.align.xg_file_name = '1.xg';
% config.stack.align.margin = 200;
% config.stack.align.roi = [];
% (f) Stack resolution parameters for area and boundary histograms. These
% are to be used if the models having been trained on images of resolution
% different from the current one. This might improve the results,
% especially during boot-strapping.
config.stack.length_factor = 1;
config.stack.area_factor = config.stack.length_factor.^2;
% () Set the scale at which segmentation should be performed. This is kept
% constant for all superpixel and superpixel-to-seg stages. When stitching
% the segmentations, the label maps can be further downsampled if needed,
% but not artificially upsampled. During linkage and proofreading, the
% segmentations are upsampled to original resolution.
% config.stack.segmentation_scale = 2;
% config.stack.segmentation_scaling_method = 'bilinear';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) name for the reconstruction
config.reconstruction.name = 'analysis1';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Mitochondria
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Training image name:
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
config.mitochondria.train.image_prefix = 'mitochondria.%03d';
config.mitochondria.train.image_suffix = '.tif';
% (b) image ids: : If the images are a006.tif, ..., a010.tif then the
% case_ids are 6,...,10.
config.mitochondria.train.case_ids = 1:3;
% (c) For constructing the feature vector for mitochondria detection. In
% case intensity histograms, the patch size is (2xwindows_size+1) X (2xwindows_size+1).
% The intensity histogram's bins are specified in intensity_bins
config.mitochondria.feature.type = 'heq_intensity_hist';
config.mitochondria.feature.window_sizes = [4 10 20 30]; % 15 25 35
config.mitochondria.feature.intensity_bins = 0:0.1:1;
% (d) Classifier type
config.mitochondria.model.type = 'boost';
% (e) Number of iterations of boosting
config.mitochondria.model.n_iteration = 30; % 30
% (f) Tree depth
config.mitochondria.model.tree_depth = 2; % 1

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Vesicle
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.vesicle.apply.dir = '';
config.vesicle.apply.save_suffix = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Superpixel - multi stage
%%%%%%%%%%%%%%%%%%%%%%%%%%
if(length(config.superpixel)==1)
  config.superpixel = repmat(config.superpixel, [3 1]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.1 Superpixel - agglomerative mean boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel(1).method = 'grayscale_AglMeanB';
% () thresholds on the boundary field
config.superpixel(1).f_threshold_seq = 0.01:0.01:0.05;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel(1).save_f_thresholds = 0.05;
% Whether to save the segmentation overlaps as TIF files for viewing.
config.superpixel(1).save_segmentation_overlay = true;
% () minimum boundary length (pixel units) for two segments to be merged
% should increase with increasing f_thresholds. E.g., 0 for f_thresholds
% <=0.65 and 20 otherwise
config.superpixel(1).length_threshold_seq = ...
  25*(config.superpixel(1).f_threshold_seq>0.65);
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel(1).use_mitochondria = true;
% () The confidence threshold to be applied on the mitochondria
config.superpixel(1).mitochondria_confidence_threshold = 0; 
% () The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel(1).mitochondria_erosion = 5; % 20
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
% Image filtering for boundary detection to be used for watershed
% boundary = 1 - filtered_image
config.superpixel(1).watershed_filter_version = 'db_nwo500_LDA1_mf3_ps1_neg'; % 'BELc\_tbb5\_mf7
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel(1).filter_version = 'db_nwo500_LDA1_mf3_ps1_neg';
% () whether to print messages
config.superpixel(1).is_verbose = is_verbose;
% () whether to display intermediate results
config.superpixel(1).is_verbose_figures = is_verbose_figures;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.3 Superpixel - grayscale_AglMedianB
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel(2).method = 'grayscale_AglMedianB';
% () Superpixel method
config.superpixel(2).superpixel_method = 'grayscale_AglMeanB';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
s = {};
for f = config.superpixel(1).save_f_thresholds %0.1
  s{end+1} = ...
    [ '.gs_amb_T', num2str(f,'%g'), '_', num2str(f-0.01,'%g'), ...
      '_b', num2str(25*(f>0.65),'%d'), ...
      '.db_nwo300_LDA3_mf3_ps1_0_neg_m0_d5_a40', ...
      '.db_nwo300_LDA3_mf3_ps1_0_neg']; %#ok<AGROW>
end
config.superpixel(2).superpixel_suffix = s; 
% () thresholds on the boundary field
config.superpixel(2).f_threshold_seq = 0.01:0.01:0.6;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel(2).save_f_thresholds = 0.3:0.05:0.6; % 0.45:0.01:0.5; 
% Whether to save the segmentation overlaps as TIF files for viewing.
config.superpixel(2).save_segmentation_overlay = true;
% () minimum boundary length (pixel units) for two segments to be merged
% should increase with increasing f_thresholds. E.g., 0 for f_thresholds
% <=0.65 and 20 otherwise
config.superpixel(2).length_threshold_seq = ...
  5*(config.superpixel(2).f_threshold_seq>-1);
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel(2).use_mitochondria = false;
% () The confidence threshold to be applied on the mitochondria
config.superpixel(2).mitochondria_confidence_threshold = 0; 
% () The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel(2).mitochondria_erosion = 5; % 20
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
config.superpixel(2).filter_version = 'db_nwo300_LDA3_mf3_ps1\_0_neg';
% () whether to print messages
config.superpixel(2).is_verbose = is_verbose;
% () whether to display intermediate results
config.superpixel(2).is_verbose_figures = is_verbose_figures;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.3 Superpixel - grayscale_ladder
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel(3).method = 'grayscale_ladder';
% () Superpixel method
config.superpixel(3).superpixel_method = 'grayscale_AglMedianB';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
s = {};
for f = config.superpixel(2).save_f_thresholds %0.1
  s{end+1} = ...
    [ '.gs_amdb_sp_T', num2str(f,'%g'), '_', num2str(f-0.01,'%g'), ...
      '_b', num2str(25*(f>0.65),'%d'), ...
      '.gs_amb_T0.05_0.04_b0', ...
      '.db_nwo300_LDA3_mf3_ps1_0_neg_m0_d5_a40', ...
      '.db_nwo300_LDA3_mf3_ps1_0_neg', ...
      '.db_nwo300_LDA3_mf3_ps1_0_neg']; %#ok<AGROW>
end
config.superpixel(3).superpixel_suffix = s; 
% () thresholds on the boundary fieldss
config.superpixel(3).f_thresholds = 0.00001;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel(3).save_f_thresholds = config.superpixel(3).f_thresholds;
% Whether to save the segmentation overlaps as TIF files for viewing.
config.superpixel(3).save_segmentation_overlay = true;
% () minimum area of a segment - below this they merged with neighbors
config.superpixel(3).area_thresholds = 300;
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel(3).use_mitochondria = false;
% () The confidence threshold to be applied on the mitochondria
config.superpixel(3).mitochondria_confidence_threshold = 0.5; 
% () The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel(3).mitochondria_erosion = 2; % 20
% * Minimum of a thresholded connected component for it to be considered
% mitochondria
config.superpixel(3).mitochondria_min_area = 40;
% () Whether to use vesicle detection to suppress false boundaries
config.superpixel(3).use_vesicle = false;
% () The threshold for the vesicle detection confidence
config.superpixel(3).vesicle_threshold = 0;
% () Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.
config.superpixel(3).vesicle_dilation = 0;
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel(3).filter_version = 'db_nwo300_LDA3_mf3_ps1\_0_neg';
% () whether to print messages
config.superpixel(3).is_verbose = is_verbose;
% () whether to display intermediate results
config.superpixel(3).is_verbose_figures = is_verbose_figures;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 14. Evaluate segment boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%
% () Ground truth segmentation to be used for evaluation.
% Format similar to that for volume reconstruction evaluation.
% seg is a 3D matrix of body labels for the voxels.
% config.eval_segment_boundary.groundtruth_file = ...
%   'region.fs.12.16.fs.12_13_14_15_16/analysis1/satoko.031009/proofread.a012.cropped001.seg.mat';
config.eval_segment_boundary.proofread_name = 'satoko.031009';
% () Lsit of segments in the ground-truth that may be over-segmented without
% penalty. Use this for irrelevant structurs , e.g., glia.
config.eval_segment_boundary.oversegment_ignore_groundtruth_label = [];
% () Method to be evaluated
% and Parameters for the method to be evaluated in the segmentation suffix
% Specify a list of parameters. Provide the ANSI C
% compatible printf format string and the list of ids to be printed with
% the format string.
% config.eval_segment_boundary.method = 'grayscale_ladder';
% config.eval_segment_boundary.seg_suffix_format = ...
% [ '.gs_l_sp_T0.05_L300', ...
%   '.gs_amdb_sp%s', ...
%   '.gs_amb_T0.05_0.04_b0', ...
%   '.db_nwo500_LDA1_mf3_ps1_neg_m0_d5_a40', ...
%   '.db_nwo500_LDA1_mf3_ps1_neg', ...
%   '.db_nwo500_LDA1_mf3_ps1_neg', ...
%   '.db_nwo500_LDA1_mf3_ps1_neg'];
% config.eval_segment_boundary.seg_suffix_id = {''};
% Parameters used for reconstructing the ministack - for baseline.
% () Maximum matching distance between ground-truth and automatic
% segmentation boundaries. Depends upon resolution. Keeping it very small
% results in zig-zag boundaries being penalized even though the overall
% segmentation is OK.
config.eval_segment_boundary.max_match_dist = 15;
% () Verbose
config.eval_segment_boundary.is_verbose = is_verbose;
% () Verbose
config.eval_segment_boundary.is_verbose_figures = is_verbose_figures;
% % () Mode of evaluation: 'pr-curve', 'match-save'
% % pr-curve: A precision-recall is generated, typically for a range of
% % segmenation parameters specified using seg_suffix_id.
% % match-save: Saves the bipartite matching results to .tif images and
% % creates a .tex file. Compiling this file using pdflatex produces a pdf of
% % the matching results. This enables detailed inspection of the evaluation.
% % config.eval_segment_boundary.mode = 'pr-curve';
% config.eval_segment_boundary.mode = 'match-save';
% type of marker for the pr-curve
config.eval_segment_boundary.marker = 'd-';

return
end
