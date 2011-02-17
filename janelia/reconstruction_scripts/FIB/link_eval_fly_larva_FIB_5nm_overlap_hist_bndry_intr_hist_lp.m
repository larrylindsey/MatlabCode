function link_eval_fly_larva_FIB_5nm_overlap_hist_bndry_intr_hist_lp(module_id, case_ids)
% reconstruction_fly_larva_FIB_5x5x5nm_ssTEM_pipeline(module_id)
% Evaluate linkage in fly larva FIB 5x5x5nm stack
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
% v1  08152008  lamina J recon.
% v2  08252008  test SIFT-based alignment
% v3  10142008  modified for medulla
% v4  12112008  modified to test tile pair alignment under folds using
%                 deformable mesh
% v5  12152008  Example of reconstruction script to be used as template
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. Constants - Not to be changed by users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config = get_basic_config();

%%% Stack root directory redefinition:
% It is recommended that data be stored in
% chklovskiilab/electron_microscopy_data/. Atypically, the medulla Leginon
% stack is stored in a different directory. This was done because of the
% very large size of the data and the duration of the project.
config.stack.dir = '/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/data_em_images/';
%%% Reconstruction root directory redefinition:
% It is recommended that all intermediate files generated during
% reconstruction be stored in the user's home directory
% em_reconstruction/reconstructions/. However, the directory can redefined
% to a common share if multiple users wish to read and write on the same
% files and to avoid duplicate storage for large projects. For instance,
% uncomment the following line to save intermediate reconstruction files in
% the chklovskii/medulla/ directory:
% config.reconstruction.root_dir = '/groups/chklovskii/medulla/reconstructions/';

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
config.stack.name = 'fly_larva.FIB.5x5x5nm';   
% (b) image name:
% If the files are image_0001.tif, image_0002.tif, etc., then
% the image_prefix would be "image_%04d" and the image_suffix would be
% ".tif".
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
config.stack.image_prefix = 'dn8.%03d';  
config.stack.image_suffix = '.tif';
% (c) slice ids: If the images are a001.tif, ..., a010.tif then the
% case_ids are 1,2,...,10.
if(nargin>1)
  config.stack.case_ids = case_ids;
else
  config.stack.case_ids = 1:2;
end
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

%%%%%%
% 1.1 Precompute data-structures for subsequent routines. E.g., fold masks.
%%%%%%
%%% 1.1.fold
% Whether to compute and store fold masks for future use
config.precompute.fold.is_enabled = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Alignment parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Currently (01/08/2009) alignment routines require trakEM stye XML
% structure to specify the structure of the stack.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% name for the reconstruction
config.reconstruction.name = '09072009';

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
config.mitochondria.train.case_ids = 6:10;
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
% (a) Algorithm to be used - must match with the function called
config.superpixel.method = 'grayscale_ladder';
% (b) thresholds on the boundary field
config.superpixel.f_thresholds = 0.36;
% (c) minimum area of a segment - below this they merged with neighbors
config.superpixel.area_thresholds = 200;
% (d) Whether mitochondria are to be used to suppress false boundaries
config.superpixel.use_mitochondria = false;
% (e) The confidence threshold to be applied on the mitochondria
config.superpixel.mitochondria_confidence_threshold = 0.0;
% (f) The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel.mitochondria_erosion = -1;
% * Minimum of a thresholded connected component for it to be considered
% mitochondria
config.superpixel.mitochondria_min_area = 500;
% (g) Whether to use vesicle detection to suppress false boundaries
config.superpixel.use_vesicle = false;
% (h) The threshold for the vesicle detection confidence
config.superpixel.vesicle_threshold = 0;
% (i) Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.
config.superpixel.vesicle_dilation = 0;
% (j) Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel.filter_version = 'mf3';
config.superpixel.is_verbose = true;
config.superpixel.is_verbose_figures = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7 Choose a superpixel segmentation parameter for further processing
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two options:
% --- (a) One parameter for all images ---
% Superpixel method
config.superpixel_choose.choice.method = 'grayscale_ladder';
% use superpixel segmentation parameters from previous step
% E.g., for grayscale-ladder it would be of the form
% .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria m<threshold_d<amount of erosion>>  
config.superpixel_choose.choice.seg_suffix = ...
  '.gs_l_T0.36_L200.mf3.mf3';
% --- (b) For each image choose a parameter from a set ---
% % Use a GUI to choose between a set of parameters for each image.
% config.superpixel_choose.choice.method = '#(SP_CHOICE)';
% config.superpixel_choose.choice.seg_suffix = '#(SP_CHOICE)';
% % * This is a list of superpixel parameters. Each element is a struct with
% % fields:
% %   .method       the segmentation method, e.g., 'grayscale_ladder'.
% %   .seg_suffix   parameters for segmentation, similar to those used in
% %                   superpixel_suffix in reconstruction scripts.
% config.superpixel_choose.param(1).method = 'grayscale_AglBIFL';
% config.superpixel_choose.param(1).seg_suffix = ...
%   '.gs_l_T0.44_L600.v_heq_mf10_m0_d-1';
% config.superpixel_choose.param(2).method = 'grayscale_AglBIFL';
% config.superpixel_choose.param(2).seg_suffix = ...
%   '.gs_l_T0.52_L600.v_heq_mf10_m0_d-1';
% config.superpixel_choose.param(3).method = 'grayscale_AglBIFL';
% config.superpixel_choose.param(3).seg_suffix = ...
%   '.gs_l_T0.56_L600.v_heq_mf10_m0_d-1';
% config.superpixel_choose.param(4).method = 'grayscale_AglBIFL';
% config.superpixel_choose.param(4).seg_suffix = ...
%   '.gs_l_T0.64_L600.v_heq_mf10_m0_d-1';
% % * Desired decision tree: This specifies the manner in which choices are
% % to be presented to the user. Two options:
% % (1) use one of the predefined decision trees (see get_basic_config.m), OR
% % (2) define a custom decision tree: See get_basic_config, Sec. 12 for examples.
% config.superpixel_choose.decision_tree = config.decision_tree_predef_binary_4_choice;
% % Whether to display the image after histogram equalization
% config.superpixel_choose.is_image_histeq = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8. Superpixel to segment
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel_2_seg.method = 'grayscale_AglMeanB';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
config.superpixel_2_seg.superpixel_method = ''; 
config.superpixel_2_seg.superpixel_suffix = ''; 
% (c) thresholds on the boundary field
config.superpixel_2_seg.f_threshold_seq = 0.36:0.01:0.55;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel_2_seg.save_f_thresholds = 0.55:0.01:0.62;
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
config.superpixel_2_seg.filter_version = 'mf3';
% () whether to display intermediate results
config.superpixel_2_seg.is_verbose = true;
config.superpixel_2_seg.is_verbose_figures = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9 Choose a segmentation parameter for further processing
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two options:
% --- (a) One parameter for all images ---
% Superpixel method
config.segmentation_choose.choice.method = 'grayscale_AglMeanB';
% use segmentation parameters from previous step
% E.g., for grayscale-ladder it would be of the form
% .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria m<threshold_d<amount of erosion>>  
config.segmentation_choose.choice.seg_suffix = ...
  [ '.gs_amb_sp_T0.55_0.54_b0', ...
    config.superpixel_choose.choice.seg_suffix, ...
    '.mf3'];
% --- (b) For each image choose a parameter from a set ---
% % Use a GUI to choose between a set of parameters for each image.
% config.superpixel_choose.choice.method = '#(SEG_CHOICE)';
% config.superpixel_choose.choice.seg_suffix = '#(SEG_CHOICE)';
% % * This is a list of superpixel parameters. Each element is a struct with
% % fields:
% %   .method       the segmentation method, e.g., 'grayscale_ladder'.
% %   .seg_suffix   parameters for segmentation, similar to those used in
% %                   superpixel_suffix in reconstruction scripts.
% config.superpixel_choose.param(1).method = 'grayscale_ladder';
% config.superpixel_choose.param(1).seg_suffix = ...
%   '.gs_l_sp_T0.58_L600.gs_l_T0.56_L600.v_heq_mf10_m0_d-1.v_heq_mf10';
% config.superpixel_choose.param(2).method = 'grayscale_ladder';
% config.superpixel_choose.param(2).seg_suffix = ...
%   '.gs_l_sp_T0.58_L600.gs_l_T0.56_L600.v_heq_mf10_m0_d-1.v_heq_mf10';
% config.superpixel_choose.param(3).method = 'grayscale_ladder';
% config.superpixel_choose.param(3).seg_suffix = ...
%   '.gs_l_sp_T0.58_L600.gs_l_T0.56_L600.v_heq_mf10_m0_d-1.v_heq_mf10';
% config.superpixel_choose.param(4).method = 'grayscale_ladder';
% config.superpixel_choose.param(4).seg_suffix = ...
%   '.gs_l_sp_T0.58_L600.gs_l_T0.56_L600.v_heq_mf10_m0_d-1.v_heq_mf10';
% % * Desired decision tree: This specifies the manner in which choices are
% % to be presented to the user. Two options:
% % (1) use one of the predefined decision trees (see get_basic_config.m), OR
% % (2) define a custom decision tree: See get_basic_config, Sec. 12 for examples.
% config.superpixel_choose.decision_tree = config.decision_tree_predef_binary_4_choice;
% % Whether to display the image after histogram equalization
% config.superpixel_choose.is_image_histeq = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Align segmentations of tiles within a section and correspond
% overlapping segments
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11. 3D linkage graph - training and generation
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D annotation file to be used for training
config.linkage.train.manual_annotation_file = ...
  '~/research/em_reconstruction_pipeline/data_em_images/fly_larva.FIB.5x5x5nm/manual_annotations/proofread.victor.wcat.mat';
% () Stack resolution parameters for area and boundary histograms. These
% are to be used if the models having been trained on images of resolution
% different from the current one. This might improve the results,
% especially during boot-strapping.
config.stack.length_factor = 1;
config.stack.area_factor = config.stack.length_factor.^2;
% (e) Type of classifier
config.linkage.model.type = 'boost_lp_p3';
% Boost classifier for overlap area (across sections)
config.linkage.model.overlap_hist.n_iteration = 30;
config.linkage.model.overlap_hist.tree_depth = 2;
config.linkage.model.overlap_hist.merge_error_relative_weight = 20;
% Boost classifier for overlap area (across sections)
config.linkage.model.boundary_hist.n_iteration = 40;
config.linkage.model.boundary_hist.tree_depth = 2;
config.linkage.model.boundary_hist.merge_error_relative_weight = 20;
% (b) Version name for the linkage model
config.linkage.train.save_suffix = [ ...
  num2str(config.linkage.model.overlap_hist.n_iteration), '.', ...
  num2str(config.linkage.model.overlap_hist.tree_depth), '.', ...
  num2str(config.linkage.model.overlap_hist.merge_error_relative_weight), '.', ...
  num2str(config.linkage.model.boundary_hist.n_iteration), '.', ...
  num2str(config.linkage.model.boundary_hist.tree_depth), '.', ...
  num2str(config.linkage.model.boundary_hist.merge_error_relative_weight)];
% (c) Type for feature to be used for linkage
config.linkage.feature.type = 'overlap_hist_boundary_interior_hist';
config.linkage.feature.boundary_hist.bin_size = 5;
config.linkage.feature.interior_hist.bin_size = 10;
% (h) The model version to be used when applying the linkage on the stack
config.linkage.apply.model_suffix = config.linkage.train.save_suffix;

%%%
% For evaluating the linkage module
%%%
% % The data is evaluated with a round-robin protocol. It is split into 2
% % sets. In each round, one set is used for testing and the rest are used
% % for training.
% Segmentation maps to be used for learning and testing. E.g., the
% ground-truth segmentation maps may be modified using Level-sets.
config.linkage.evaluate.segmentation_file = ...
  '/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/fly_larva.FIB.5x5x5nm/2D_segmentation_results/segmentation_for_linkage_evaluation/seg.wcat.gs_amb_sp_T0.55_0.54_b0.gs_l_T0.36_L200.mf3.mf3.mf3.mat';
% Whether to train on the segmentation maps in segmentation_file or on the
% manual_annotation.
config.linkage.evaluate.is_trained_on_manual_annotation = false;
% config.linkage.evaluate.k_fold = 2;
% Linkage thresholds to used for generating the
% false-acceptance/false-rejection curves.
config.linkage.evaluate.link_thresholds = 0:0.1:2;
% whether to display the link decisions
config.linkage.evaluate.display_links = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12. Dump to proofreader
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Which proofreader is being used. E.g., matlab_gui, Raveler.
config.proofreader.method = 'Raveler';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. Reconstruction Routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pipeline(module_id, config);

return

end