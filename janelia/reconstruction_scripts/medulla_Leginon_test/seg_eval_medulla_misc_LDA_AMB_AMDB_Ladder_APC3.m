function config = seg_eval_medulla_misc_LDA_AMB_AMDB_Ladder_APC3(config, case_ids, is_verbose, is_verbose_figures)
% seg_eval_medulla_misc_LDA_ladder(config, case_ids, is_verbose,
% is_verbose_figures)
% Evaluate LDA on medulla misc. dataset.
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
% () stack name
config.stack.name = 'medulla.HPF.Leginon.3500x.zhiyuan.fall2008';   
% () image name:
% Option 1: If the files are image_0001.tif, image_0002.tif, etc., then
% the image_prefix would be "image_%04d" and the image_suffix would be
% ".tif".
% For files named as a001.tif, a002.tif, etc., the image_prefix would be
% "a%03d" and the image_suffix would be ".tif". 
% If the images are a001.tif, ..., a010.tif then the
% case_ids would be 1:10.
config.stack.image_prefix = 'a.%03d';
config.stack.image_suffix = '.tif';
% Option 2: If the images correspond to a trakEM project then
% * specify the trakEM file structure xml file name
% config.stack.image_structure = 'eval_data_596_605_1tile.xml';
% config.stack.image_structure = 'crop4_global_alignment_0161_0860.xml';
% The z-sections within the stack that form the region. If unspecified,
% this defaults to the config.stack.case_ids. Together with
% config.stack.image_structure, config.region.case_ids decides the region
% directory and hence the data_to_be_proofread directory. This should be
% kept constant a particular region's reconstruction but MUST BE CHANGED
% when shifting to a different region.
if(nargin>1 && ~isempty(case_ids))
  config.stack.case_ids = case_ids;
else
  config.stack.case_ids = 12:16;
end
% () ROI: Region Of Interest for the stack. If the entire slice image is
% to be processed then assign empty ([]). Otherwise get_roi(config);
config.stack.roi = [];
% () Stack resolution parameters for area and boundary histograms. These
% are to be used if the models having been trained on images of resolution
% different from the current one. This might improve the results,
% especially during boot-strapping.
config.stack.length_factor = 1;
config.stack.area_factor = config.stack.length_factor.^2;
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
config.superpixel(1).watershed_filter_version = 'db_nwo700_cbf2_LDA7_mf7_ps1\_0.25_abf_neg';
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel(1).filter_version = 'db_nwo700_cbf2_LDA7_mf7_ps1\_0.25_abf_neg';
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
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg_m0_d5_a40', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg']; %#ok<AGROW>
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
config.superpixel(2).filter_version = 'db_nwo700_cbf2_LDA7_mf7_ps1\_0.25_abf_neg';
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
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg_m0_d5_a40', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg']; %#ok<AGROW>
end
config.superpixel(3).superpixel_suffix = s; 
% () thresholds on the boundary fieldss
config.superpixel(3).f_thresholds = 0.3;
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
config.superpixel(3).filter_version = 'db_nwo700_cbf2_LDA7_mf7_ps1\_0.25_abf_neg';
% () whether to print messages
config.superpixel(3).is_verbose = is_verbose;
% () whether to display intermediate results
config.superpixel(3).is_verbose_figures = is_verbose_figures;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.1 Superpixel - prune-classify
%%%%%%%%%%%%%%%%%%%%%%%%%%
if(length(config.superpixel_2_seg)==1)
  config.superpixel_2_seg = repmat(config.superpixel_2_seg, [3 1]);
end
% (a) Algorithm to be used - must match with the function called
config.superpixel_2_seg(1).method = 'prune_classify';
% () Superpixel method
config.superpixel_2_seg(1).superpixel_method = 'grayscale_ladder';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
s = {};
for f = config.superpixel(2).save_f_thresholds %0.1
  s{end+1} = ...
    [ '.gs_l_sp_T0.3_L300', ...
      '.gs_amdb_sp_T', num2str(f,'%g'), '_', num2str(f-0.01,'%g'), ...
      '_b', num2str(5*(f>-1),'%d'), ...
      '.gs_amb_T0.05_0.04_b0', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg_m0_d5_a40', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg']; %#ok<AGROW>
end
config.superpixel_2_seg(1).superpixel_suffix = s; 
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
config.superpixel_2_seg(1).prune_classify.merge_criterion_param.patch_size = 3;
% For region_max: the threshold on the maximum value observed within a
% patch
config.superpixel_2_seg(1).prune_classify.merge_criterion_param.min_threshold = 0.35;
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
config.superpixel_2_seg(1).filter_version = 'db_nwo700_cbf2_LDA7_mf7_ps1\_0.25_abf_neg';
% () whether to print messages
config.superpixel_2_seg(1).is_verbose = is_verbose;
% () whether to display intermediate results
config.superpixel_2_seg(1).is_verbose_figures = is_verbose_figures;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.2 Superpixel - prune-classify
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel_2_seg(2).method = 'prune_classify';
% () Superpixel method
config.superpixel_2_seg(2).superpixel_method = 'prune_classify';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
s = {};
for f = config.superpixel(2).save_f_thresholds(3) %0.1
  s{end+1} = ...
    [ '.apc_sp_region_max_3_0.35', ...
      '.gs_l_sp_T0.3_L300', ...
      '.gs_amdb_sp_T', num2str(f,'%g'), '_', num2str(f-0.01,'%g'), ...
      '_b', num2str(5*(f>-1),'%d'), ...
      '.gs_amb_T0.05_0.04_b0', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg_m0_d5_a40', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg']; %#ok<AGROW>
end
config.superpixel_2_seg(2).superpixel_suffix = s; 
% () criterion used for merging.
% Options:
% - 'min': minimum value encountered along a boundary
% - 'region_max': take patches along the boundary and compute their
%   maximum. Merge is one of these is below a threshold. Looks for large
%   gaps in the boundary confidence.
% - 'boundary_hist_boost': compute boundary and interior histograms to
%   represent pairs adjacent segments and classify using boosting.
config.superpixel_2_seg(2).prune_classify.merge_criterion = 'boundary_hist_wekarf';
% For boundary_hist_boost: the threshold on the boost confidence for
% merging.s
config.superpixel_2_seg(2).prune_classify.merge_criterion_param.threshold = 0.7;
% Whether to save the segmentation overlaps as TIF files for viewing.
config.superpixel_2_seg(2).save_segmentation_overlay = true;
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel_2_seg(2).use_mitochondria = false;
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
config.superpixel_2_seg(2).filter_version = 'db_nwo700_cbf2_LDA7_mf7_ps1\_0.25_abf_neg';
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
s = {};
for f = config.superpixel(2).save_f_thresholds %0.1
  s{end+1} = ...
    [ '.apc_sp_region_max_3_0.35', ...
      '.gs_l_sp_T0.3_L300', ...
      '.gs_amdb_sp_T', num2str(f,'%g'), '_', num2str(f-0.01,'%g'), ...
      '_b', num2str(5*(f>-1),'%d'), ...
      '.gs_amb_T0.05_0.04_b0', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg_m0_d5_a40', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg']; %#ok<AGROW>
end
config.superpixel_2_seg(3).superpixel_suffix = s; 
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
config.superpixel_2_seg(3).filter_version = 'db_nwo700_cbf2_LDA7_mf7_ps1\_0.25_abf_neg';
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
%     '.BELc_5_mf7_e2_d4_f', ...
%     '.BELc_5_mf7_e2_d4_f'];
config.segmentation_choose.choice.method = 'prune_classify';
config.segmentation_choose.choice.seg_suffix = ...
  [ '.apc_sp_region_max_3_0.35', ...
    '.gs_l_sp_T0.3_L300', ...
    '.gs_amdb_sp_T0.5_0.49_b5', ...
    '.gs_amb_T0.05_0.04_b0', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg_m0_d5_a40', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg'];
% --- (b) For each image choose a parameter from a set ---
% % Use a GUI to choose between a set of parameters for each image.
% --- (c) Set a parameter for a batch of sections ---
% % Useful if blocks of sections have different parameters.
% config.segmentation_choose.set_choice.method = 'prune_classify';
% config.segmentation_choose.set_choice.seg_suffix = ...
%   [ '.apc_sp_region_max_2_0.2', ...
%     config.superpixel_choose.choice.seg_suffix, ...
%     '.db_nwo650_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg'];

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
config.linkage.feature.filter_version = 'db_nwo700_cbf2_LDA7_mf7_ps1\_0.25_abf_neg';
% Feature suffix for the file name
config.linkage.feature.suffix = ['.', num2str(config.linkage.feature.bin_size), ...
  '.', config.linkage.feature.filter_version];
% Use a modified version of segmentation for the linkage - specify prefix
% and suffix to the segmentation parameter string for retrieving the files
config.linkage.modified_segment.method = 'shrink_levelset';
config.linkage.modified_segment.prefix = '.skl_sp_mul5_i25';
config.linkage.modified_segment.suffix = '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg_m0.5_d2_a40_f';
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
config.eval_segment_boundary.seg_suffix_format = ...
    [ '.apc_sp_boundary_hist_wekarf_0.7', ...
      '.apc_sp_region_max_3_0.35', ...
      '.gs_l_sp_T0.3_L300', ...
      '.gs_amdb_sp_T%s', ...
      '.gs_amb_T0.05_0.04_b0', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg_m0_d5_a40', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg'];
s = {};
for f = config.superpixel(2).save_f_thresholds(3)
  s{end+1} = [num2str(f,'%g'), '_', num2str(f-0.01,'%g'), ...
      '_b', num2str(5*(f>-1),'%d')]; %#ok<AGROW>
end
config.eval_segment_boundary.seg_suffix_id = s;
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
% type of marker for the pr-curve
config.eval_segment_boundary.marker = 'd-';

return
end
