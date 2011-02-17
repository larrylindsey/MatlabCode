function config = seg_eval_medulla_misc_LDA_apc_bndry_hist_boost(config, case_ids, is_verbose, is_verbose_figures)
% linkage_eval_medulla_596_605_s_LDA_apc_bndry_hist_boost(module_id)
% Evaluate linkage on medulla.HPF.Leginon.3500x.zhiyuan.fall2008.
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. Constants - Not to be changed by users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Whether to use hashing when generating file names during storage and
% % retrieval. This reduces the filename lengths.
% global config_global
% config_global.storage.file_name_is_hashed__ = false;

config.stack.dir = '/groups/chklovskii/medulla/';

% For temporarily saving the superpixel-to-segment and segment-to-body
% maps.
config.to_be_proofread.root_dir = '~/research/em_reconstruction_pipeline/data_to_be_proofread/';

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
% Option 1: If the files are image_0001.tif, image_0002.tif, etc., then
% the image_prefix would be "image_%04d" and the image_suffix would be
% ".tif".
% For files named as a001.tif, a002.tif, etc., the image_prefix would be
% "a%03d" and the image_suffix would be ".tif". 
% If the images are a001.tif, ..., a010.tif then the
% case_ids would be 1:10.
config.stack.image_prefix = 'a.%03d';
config.stack.image_suffix = '.tif';
if(nargin>1 && ~isempty(case_ids))
  config.stack.case_ids = case_ids;
else
  config.stack.case_ids = 12:16;
end
% % Option 2: If the images correspond to a trakEM project then
% % * specify the trakEM file structure xml file name
% config.stack.image_structure = 'eval_data_596_605_1tile.xml';
% % The z-sections within the stack that form the region. If unspecified,
% % this defaults to the config.stack.case_ids. Together with
% % config.stack.image_structure, config.region.case_ids decides the region
% % directory and hence the data_to_be_proofread directory. This should be
% % kept constant a particular region's reconstruction but MUST BE CHANGED
% % when shifting to a different region.
% config.region.case_ids = 599:605;
% % * slice ids: specify the z-plane numbers.
% if(nargin>1 && ~isempty(case_ids))
%   config.stack.case_ids = case_ids;
% else
%   config.stack.case_ids = 596:600;
% end
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
% config.stack.align.is_prealigned = true;
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
% () Fold detection and consideration during reconstruction
% Whether folds are considered during segmentation
config.stack.fold.is_considered_in_segmentation = false;
% Whether folds are considered during alignment
config.stack.fold.is_considered_in_alignment = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) name for the reconstruction
config.reconstruction.name = 'link_eval_ls_lp';

%%%%%%
% 1.1 Precompute data-structures for subsequent routines. E.g., fold masks.
%%%%%%
%%% 1.1.fold
% Whether to compute and store fold masks for future use
config.precompute.fold.is_enabled = true;
% Whether to save tif files of patches
config.precompute.fold.save_patch_tifs = true;
% Minimum area of a component in the thresholded image for it to be
% considered as a fold.
config.precompute.fold.min_fold_area = 1000;
% Suffix for storing fold masks
config.fold.save_suffix = '.foldmask';

%%% 1.1.filtered_image
% Precompute filtered images
config.precompute.filtered_image.filter_version = 'db_nwo650_cbf2_LDA7_mf7_ps1.5\_0.2_abf_mito0\_d5\_e5\_a400_neg';


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
config.superpixel(1).watershed_filter_version = 'db_nwo650_cbf2_LDA7_mf7_ps1.5\_0.2_abf_neg';
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel(1).filter_version = 'db_nwo650_cbf2_LDA7_mf7_ps1.5\_0.2_abf_neg';
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
      '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg_m0_d5_a40', ...
      '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg']; %#ok<AGROW>
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
config.superpixel(2).filter_version = 'db_nwo650_cbf2_LDA7_mf7_ps1.5\_0.2_abf_neg';
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
      '_b', num2str(5*(f>-1),'%d'), ...
      '.gs_amb_T0.05_0.04_b0', ...
      '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg_m0_d5_a40', ...
      '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
      '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg']; %#ok<AGROW>
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
config.superpixel(3).filter_version = 'db_nwo650_cbf2_LDA7_mf7_ps1.5\_0.2_abf_neg';
% () whether to print messages
config.superpixel(3).is_verbose = is_verbose;
% () whether to display intermediate results
config.superpixel(3).is_verbose_figures = is_verbose_figures;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7 Choose a superpixel segmentation parameter for further processing
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two options:
% --- (a) One parameter for all images ---
% % Superpixel method
% config.superpixel_choose.choice.method = 'grayscale_ladder';
% % use superpixel segmentation parameters from previous step
% % E.g., for grayscale-ladder it would be of the form
% % .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria m<threshold_d<amount of erosion>>  
% config.superpixel_choose.choice.seg_suffix = ...
%     [ '.gs_l_sp_T1e-05_L300', ...
%       '.gs_amdb_sp_T0.4_0.39_b5', ...
%       '.gs_amb_T0.05_0.04_b0', ...
%       '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg_m0_d5_a40', ...
%       '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
%       '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
%       '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg'];
% --- (b) For each image choose a parameter from a set ---
% % Use a GUI to choose between a set of parameters for each image.
% config.superpixel_choose.choice.method = '#(SP_CHOICE)';
% config.superpixel_choose.choice.seg_suffix = '#(SP_CHOICE)';
% % % % * This is a list of superpixel parameters. Each element is a struct with
% % % % fields:
% % % %   .method       the segmentation method, e.g., 'grayscale_ladder'.
% % % %   .seg_suffix   parameters for segmentation, similar to those used in
% % % %                   superpixel_suffix in reconstruction scripts.
% % % config.superpixel_choose.param(1).method = 'grayscale_AglBIFL';
% % % config.superpixel_choose.param(1).seg_suffix = ...
% % %   '.gs_abif_sp_T0.18_0.17_b2.gs_l_sp_T0.2_L60.gs_amb_T0.45_0.44_b0.v_heq_o1_m0_d5_a40_f.v_heq2.v_heq2_f.v0';
% % % config.superpixel_choose.param(2).method = 'grayscale_AglBIFL';
% % % config.superpixel_choose.param(2).seg_suffix = ...
% % %   '.gs_abif_sp_T0.18_0.17_b2.gs_l_sp_T0.2_L60.gs_amb_T0.51_0.5_b0.v_heq_o1_m0_d5_a40_f.v_heq2.v_heq2_f.v0';
% % % config.superpixel_choose.param(3).method = 'grayscale_AglBIFL';
% % % config.superpixel_choose.param(3).seg_suffix = ...
% % %   '.gs_abif_sp_T0.18_0.17_b2.gs_l_sp_T0.2_L60.gs_amb_T0.55_0.54_b0.v_heq_o1_m0_d5_a40_f.v_heq2.v_heq2_f.v0';
% % % config.superpixel_choose.param(4).method = 'grayscale_AglBIFL';
% % % config.superpixel_choose.param(4).seg_suffix = ...
% % %   '.gs_abif_sp_T0.18_0.17_b2.gs_l_sp_T0.2_L60.gs_amb_T0.59_0.58_b0.v_heq_o1_m0_d5_a40_f.v_heq2.v_heq2_f.v0';
% % * Desired decision tree: This specifies the manner in which choices are
% % to be presented to the user. Two options:
% % (1) use one of the predefined decision trees (see get_basic_config.m), OR
% % (2) define a custom decision tree: See get_basic_config, Sec. 12 for examples.
% config.superpixel_choose.decision_tree = config.decision_tree_predef_binary_4_choice;
% config.superpixel_choose.decision_tree = config.decision_tree_predef_binary_4_choice_bias_4;
% % Whether to display the image after histogram equalization
% config.superpixel_choose.is_image_histeq = true;
% % --- (c) Set a parameter for a batch of sections ---
% % Useful if blocks of sections have different parameters.
config.superpixel_choose.choice.method = 'grayscale_ladder';
config.superpixel_choose.choice.seg_suffix = ...
  [ '.gs_l_sp_T1e-05_L300', ...
      '.gs_amdb_sp_T0.4_0.39_b5', ...
      '.gs_amb_T0.05_0.04_b0', ...
      '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg_m0_d5_a40', ...
      '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
      '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
      '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg'];

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.1 Superpixel - prune-classify
%%%%%%%%%%%%%%%%%%%%%%%%%%
if(length(config.superpixel_2_seg)==1)
  config.superpixel_2_seg = repmat(config.superpixel_2_seg, [3 1]);
end
% (a) Algorithm to be used - must match with the function called
config.superpixel_2_seg(1).method = 'prune_classify';
% () Superpixel method
config.superpixel_2_seg(1).superpixel_method = '';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
config.superpixel_2_seg(1).superpixel_suffix = ''; 
% () criterion used for merging.
% Options:
% - 'min': minimum value encountered along a boundary
% - 'region_max': take patches along the boundary and compute their
%   maximum. Merge is one of these is below a threshold. Looks for large
%   gaps in the boundary confidence.
config.superpixel_2_seg(1).prune_classify.merge_criterion = 'region_max';
% Parameters specific to the chosen merge criterion
% For region_max: the radius of the path. The patch's dimension is
% 2*patch_size+1
config.superpixel_2_seg(1).prune_classify.merge_criterion_param.patch_size = 2;
% For region_max: the threshold on the maximum value observed within a
% patch
config.superpixel_2_seg(1).prune_classify.merge_criterion_param.min_threshold = 0.2;
% Whether to save the segmentation overlaps as TIF files for viewing.
config.superpixel_2_seg(1).save_segmentation_overlay = true;
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel_2_seg(1).use_mitochondria = false;
% () The confidence threshold to be applied on the mitochondria
config.superpixel_2_seg(1).mitochondria_confidence_threshold = 0.5; 
% () The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel_2_seg(1).mitochondria_erosion = 2; % 20
% * Minimum of a thresholded connected component for it to be considered
% mitochondria
config.superpixel_2_seg(1).mitochondria_min_area = 40;
% () Whether to use vesicle detection to suppress false boundaries
config.superpixel_2_seg(1).use_vesicle = false;
% () The threshold for the vesicle detection confidence
config.superpixel_2_seg(1).vesicle_threshold = 0;
% () Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.
config.superpixel_2_seg(1).vesicle_dilation = 0;
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel_2_seg(1).filter_version = 'db_nwo650_cbf2_LDA7_mf7_ps1.5\_0.2_abf_neg';
% () whether to print messages
config.superpixel_2_seg(1).is_verbose = is_verbose;
% () whether to display intermediate results
config.superpixel_2_seg(1).is_verbose_figures = is_verbose_figures;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.1 Superpixel - prune-classify
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel_2_seg(2).method = 'prune_classify';
% () Superpixel method
config.superpixel_2_seg(2).superpixel_method = 'prune_classify';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
config.superpixel_2_seg(2).superpixel_suffix = ...
  [ '.apc_sp_region_max_2_0.2', ...
    config.superpixel_choose.choice.seg_suffix, ...
    '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg'];
% () criterion used for merging.
% Options:
% - 'min': minimum value encountered along a boundary
% - 'region_max': take patches along the boundary and compute their
%   maximum. Merge is one of these is below a threshold. Looks for large
%   gaps in the boundary confidence.
% - 'boundary_hist_boost': compute boundary and interior histograms to
%   represent pairs adjacent segments and classify using boosting.
config.superpixel_2_seg(2).prune_classify.merge_criterion = 'boundary_hist_boost';
% For boundary_hist_boost: the threshold on the boost confidence for
% merging.s
config.superpixel_2_seg(2).prune_classify.merge_criterion_param.threshold = 2.5;
% Whether to save the segmentation overlaps as TIF files for viewing.
config.superpixel_2_seg(2).save_segmentation_overlay = true;
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel_2_seg(2).use_mitochondria = true;
% () The confidence threshold to be applied on the mitochondria
config.superpixel_2_seg(2).mitochondria_confidence_threshold = 0.5; 
% () The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel_2_seg(2).mitochondria_erosion = 2; % 20
% * Minimum of a thresholded connected component for it to be considered
% mitochondria
config.superpixel_2_seg(2).mitochondria_min_area = 40;
% () Whether to use vesicle detection to suppress false boundaries
config.superpixel_2_seg(2).use_vesicle = false;
% () The threshold for the vesicle detection confidence
config.superpixel_2_seg(2).vesicle_threshold = 0;
% () Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.
config.superpixel_2_seg(2).vesicle_dilation = 0;
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel_2_seg(2).filter_version = 'db_nwo650_cbf2_LDA7_mf7_ps1.5\_0.2_abf_neg';
% () whether to print messages
config.superpixel_2_seg(2).is_verbose = is_verbose;
% () whether to display intermediate results
config.superpixel_2_seg(2).is_verbose_figures = is_verbose_figures;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.3 Superpixel - prune-classify
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel_2_seg(3).method = 'shrink_levelset';
% () Superpixel method
config.superpixel_2_seg(3).superpixel_method = 'prune_classify';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
config.superpixel_2_seg(3).superpixel_suffix = ...
  [ '.apc_sp_region_max_2_0.2', ...
    config.superpixel_choose.choice.seg_suffix, ...
    '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg']; 
% Whether to save the segmentation overlaps as TIF files for viewing.
config.superpixel_2_seg(3).save_segmentation_overlay = true;
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel_2_seg(3).use_mitochondria = true;
% () The confidence threshold to be applied on the mitochondria
config.superpixel_2_seg(3).mitochondria_confidence_threshold = 0.5; 
% () The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel_2_seg(3).mitochondria_erosion = 2; % 20
% * Minimum of a thresholded connected component for it to be considered
% mitochondria
config.superpixel_2_seg(3).mitochondria_min_area = 40;
% () Whether to use vesicle detection to suppress false boundaries
config.superpixel_2_seg(3).use_vesicle = false;
% () The threshold for the vesicle detection confidence
config.superpixel_2_seg(3).vesicle_threshold = 0;
% () Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.
config.superpixel_2_seg(3).vesicle_dilation = 0;
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel_2_seg(3).filter_version = 'db_nwo650_cbf2_LDA7_mf7_ps1.5\_0.2_abf_neg';
% () whether to print messages
config.superpixel_2_seg(3).is_verbose = is_verbose;
% () whether to display intermediate results
config.superpixel_2_seg(3).is_verbose_figures = is_verbose_figures;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9 Choose a segmentation parameter for further processing
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two options:
% --- (a) One parameter for all images ---
% Segmentation method
% use segmentation parameters from previous step
% E.g., for grayscale-ladder it would be of the form
% .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria m<threshold_d<amount of erosion>>  
% config.segmentation_choose.choice.seg_suffix = ...
%   [ '.gs_amdb_sp_T0.08_0.07_b0', ...
%     '.gs_amb_sp_T0.11_0.1_b0', ...
%     config.superpixel_choose.choice.seg_suffix, ...
%     '.BELc_5_mf7_e2_d4', ...
%     '.BELc_5_mf7_e2_d4'];
% config.segmentation_choose.choice.method = '#(SEG_CHOICE)';
% config.segmentation_choose.choice.seg_suffix = '#(SEG_CHOICE)';
% --- (b) For each image choose a parameter from a set ---
% % Use a GUI to choose between a set of parameters for each image.
% % * This is a list of superpixel parameters. Each element is a struct with
% % fields:
% %   .method       the segmentation method, e.g., 'grayscale_ladder'.
% %   .seg_suffix   parameters for segmentation, similar to those used in
% %                   superpixel_suffix in reconstruction scripts.
% config.superpixel_choose.param(1).method = 'grayscale_AglBIFL';
% config.superpixel_choose.param(1).seg_suffix = ...
%   '.gs_abif_sp_T0.18_0.17_b2.gs_l_sp_T0.2_L60.gs_amb_T0.63_0.62_b0.v_heq_o1_m0_d5_a40_f.v_heq2.v_heq2_f.v0';
% 
% config.superpixel_choose.param(2).method = 'grayscale_AglBIFL';
% config.superpixel_choose.param(2).seg_suffix = ...
%   '.gs_abif_sp_T0.16_0.15_b2.gs_l_sp_T0.2_L60.gs_amb_T0.63_0.62_b0.v_heq_o1_m0_d5_a40_f.v_heq2.v_heq2_f.v0';
% 
% config.superpixel_choose.param(3).method = 'grayscale_AglBIFL';
% config.superpixel_choose.param(3).seg_suffix = ...
%   '.gs_abif_sp_T0.14_0.13_b2.gs_l_sp_T0.2_L60.gs_amb_T0.63_0.62_b0.v_heq_o1_m0_d5_a40_f.v_heq2.v_heq2_f.v0';
% 
% config.superpixel_choose.param(4).method = 'grayscale_AglBIFL';
% config.superpixel_choose.param(4).seg_suffix = ...
%   '.gs_abif_sp_T0.12_0.11_b2.gs_l_sp_T0.2_L60.gs_amb_T0.63_0.62_b0.v_heq_o1_m0_d5_a40_f.v_heq2.v_heq2_f.v0';
% % * Desired decision tree: This specifies the manner in which choices are
% % to be presented to the user. Two options:
% % (1) use one of the predefined decision trees (see get_basic_config.m), OR
% % (2) define a custom decision tree: See get_basic_config, Sec. 12 for examples.
% config.superpixel_choose.decision_tree = config.decision_tree_predef_binary_4_choice;
% % Whether to display the image after histogram equalization
% config.superpixel_choose.is_image_histeq = true;
% --- (c) Set a parameter for a batch of sections ---
% Useful if blocks of sections have different parameters.
config.segmentation_choose.choice.method = 'prune_classify';
config.segmentation_choose.choice.seg_suffix = ...
  [ '.apc_sp_boundary_hist_boost_2.5', ...
    '.apc_sp_region_max_2_0.2', ...
    config.superpixel_choose.choice.seg_suffix, ...
    '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
    '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg'];

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Align segmentations of tiles within a section and correspond
% overlapping segments
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.align_segmentation.is_verbose = is_verbose;
config.align_segmentation.is_verbose_figures = is_verbose_figures;
config.align_segmentation.method = 'boundary_bipartite_match';
config.align_segmentation.delta_margin = 7;
config.align_segmentation.min_n_votes = 15;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11. 3D linkage graph - training and generation
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.linkage.is_verbose = is_verbose;
config.linkage.is_verbose_figures = is_verbose_figures;
% Type for feature to be used for linkage
config.linkage.feature.type = 'overlap_hist_boundary_interior_hist_shrink_LS';
% Number of bins in the intensity histogram
config.linkage.feature.bin_size = 10;
config.linkage.feature.boundary_hist.bin_size = 10;
config.linkage.feature.interior_hist.bin_size = 10;
% Filter image for linkage statistics
config.linkage.feature.filter_version = 'db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg';
% Feature suffix for the file name
config.linkage.feature.suffix = ['.', num2str(config.linkage.feature.bin_size), ...
  '.', config.linkage.feature.filter_version];
% Use a modified version of segmentation for the linkage - specify prefix
% and suffix to the segmentation parameter string for retrieving the files
config.linkage.modified_segment.method = 'shrink_levelset';
config.linkage.modified_segment.prefix = '.skl_sp_mul5_i25';
config.linkage.modified_segment.suffix = '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg_m0.5_d2_a40';
% (e) Type of classifier
config.linkage.model.type = 'boost_lp_p3';
% (f) Number of iterations of boosting
config.linkage.model.n_iteration = 40;
% (g) Tree depth
config.linkage.model.tree_depth = 2;
% (f) Number of iterations of boosting
config.linkage.model.n_iteration_insection = 40;
% (g) Tree depth
config.linkage.model.tree_depth_insection = 2;
% Model suffix for file name
config.linkage.model.suffix = [ ...
  '.', num2str(config.linkage.model.n_iteration), ...
  '.', num2str(config.linkage.model.tree_depth)];
% Ground-truth segmentation maps to be used for testing.
config.linkage.train.superpixel_method_ground_truth = 'proofread_result';
config.linkage.train.superpixel_suffix_ground_truth = '.prfrd_asp.satoko.042709';
config.linkage.train.segment_method_ground_truth = 'proofread_result';
config.linkage.train.segment_suffix_ground_truth = '.prfrd_seg.satoko.042709';
% Ground-truth linkage suffix
config.linkage.train.linkage_suffix_ground_truth = ...
  '.proofread.satoko.042709.v0#(SEG_CHOICE)';
% Suffix for linkage model in addition to feature and classifier model
config.linkage.train.save_suffix = '.satoko.042709';
% (h) The model version to be used when applying the linkage on the stack
config.linkage.apply.model_suffix = config.linkage.train.save_suffix;
% Whether to retain only subset of links, and if so the method for choosing
% the subset, options:
% (1) 'max_confidence': for each segment, retain at most one link, the one
% with maximum confidence.
% config.linkage.apply.retain_links_subset = 'max_confidence';
% Linkage threshold
config.linkage.apply.linkage_threshold = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12. Dump to proofreader
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Which proofreader is being used. E.g., matlab_gui, Raveler.
config.proofreader.method = 'Raveler'; %'matlab_gui';
% Whether to simply generate an identity ROI for the stack. If set to
% false, convex hulls are computed for all sections and a vertical overlap
% mask is computed.
config.proofreader.align_roi.is_identity = true;
% ROI to be applied before dumping to proofreader
config.proofreader.roi.xmin = 6511;
config.proofreader.roi.ymin = 7311;
config.proofreader.roi.xmax = 12510;
config.proofreader.roi.ymax = 12310;
config.proofreader.export.is_verbose_figures = is_verbose_figures;
% Annotations
% config.proofreader.annotations.method = 'proofread_result';
% config.proofreader.annotations.suffix = '.prfd_annot.concatenate_161_1460';

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
config.eval_segment_boundary.method = 'prune_classify';
% config.eval_segment_boundary.seg_suffix_format = ...
% [ '.apc_sp_region_max_2_0.2', ...
%   '.gs_l_sp_T1e-05_L300', ...
%   '.gs_amdb_sp%s', ...
%   '.gs_amb_T0.05_0.04_b0', ...
%   '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg_m0_d5_a40', ...
%   '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
%   '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
%   '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
%   '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg'];
config.eval_segment_boundary.seg_suffix_format = ...
[ '.apc_sp_boundary_hist_boost_2.5', ...
  '.apc_sp_region_max_2_0.2', ...
  '.gs_l_sp_T1e-05_L300', ...
  '.gs_amdb_sp%s', ...
  '.gs_amb_T0.05_0.04_b0', ...
  '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg_m0_d5_a40', ...
  '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
  '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
  '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
  '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
  '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg_m0.5_d2_a40'];
config.eval_segment_boundary.seg_suffix_id = {'_T0.4_0.39_b5'};
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

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 15. For evaluating the linkage module
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.linkage.evaluate.is_verbose = is_verbose;
config.linkage.evaluate.is_verbose_figures = is_verbose_figures;
% Automatic segmentation maps to be used for testing.
config.linkage.evaluate.segment_method_automatic = ...
  config.segmentation_choose.choice.method;
config.linkage.evaluate.segment_suffix_automatic = ...
  config.segmentation_choose.choice.seg_suffix;
% Automatic linkage suffix
config.linkage.evaluate.linkage_suffix_automatic = ...
  ['.boost_lp_p3', config.linkage.model.suffix, ...
  '.overlap_hist_boundary_interior_hist_shrink_LS', config.linkage.feature.suffix, ...
  '.satoko.042709', config.linkage.evaluate.segment_suffix_automatic];
% Ground-truth segmentation maps to be used for testing.
config.linkage.evaluate.superpixel_method_ground_truth = 'proofread_result';
config.linkage.evaluate.superpixel_suffix_ground_truth = '.prfrd_asp.satoko.042709';
config.linkage.evaluate.segment_method_ground_truth = 'proofread_result';
config.linkage.evaluate.segment_suffix_ground_truth = '.prfrd_seg.satoko.042709';
% Ground-truth linkage suffix
config.linkage.evaluate.linkage_suffix_ground_truth = ...
  '.proofread.satoko.042709.v0#(SEG_CHOICE)';
% Linkage thresholds to used for generating the
% false-acceptance/false-rejection curves
config.linkage.evaluate.linkage_thresholds = 0:0.1:2;

% Filter image for linkage statistics
config.split3D.filter_version = ['.', ...
  regexprep(config.precompute.filtered_image.filter_version, ...
  '(?<char>[^\\])\', '$<char>')];

return
end
