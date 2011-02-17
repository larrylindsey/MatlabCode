function link_exp_medulla_HPF_296_305_overlap_hist_boost(module_id, case_ids)
% link_exp_fly_larva_no_ua_overlap_area_boost(module_id)
% Linkage experiments: 
%
% module_id
%   type help pipeline_serial_section.m
% 
% Pratim Ghosh,
% Summer intern, Janelia Farm Research Campus, HHMI.
% Univ. of California, Santa Barbara.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  06202009  template file for linkage experiments. Example of linking
%                   by boosting on the area of overlap.
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
% config.stack.dir = '/groups/chklovskii/medulla/';
%%% Reconstruction root directory redefinition:
% It is recommended that all intermediate files generated during
% reconstruction be stored in the user's home directory
% em_reconstruction/reconstructions/. However, the directory can redefined
% to a common share if multiple users wish to read and write on the same
% files and to avoid duplicate storage for large projects. For instance,
% uncomment the following line to save intermediate reconstruction files in
% the chklovskii/medulla/ directory:
% config.reconstruction.root_dir = '/groups/chklovskii/medulla/reconstructions/';
%%% Directory for saving the to-be-proofread data for Raveler, Matlab-GUI,
%%% etc.
% config.to_be_proofread.root_dir = ...
%   '/groups/chklovskii/home/ghoshp/em_reconstruction/data_to_be_proofread/';

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
config.stack.image_prefix = 'a%03d';  
config.stack.image_suffix = '.tif';
% (c) slice ids: If the images are a001.tif, ..., a010.tif then the
% case_ids are 1,2,...,10.
if(nargin>1)
  config.stack.case_ids = case_ids;
else
  config.stack.case_ids = 3:12;
end
% (d) ROI: Region Of Interest for the stack. If the entire slice image is
% to be processed then assign empty ([]). Otherwise get_roi(config);
config.stack.roi = [];
% config.stack.roi = get_roi(config);
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
config.reconstruction.name = 'test_overlap_area';

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
config.superpixel.f_thresholds = [0.44, 0.52, 0.56, 0.64, 0.72];
% (c) minimum area of a segment - below this they merged with neighbors
config.superpixel.area_thresholds = 600;
% (d) Whether mitochondria are to be used to suppress false boundaries
config.superpixel.use_mitochondria = true;
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
config.superpixel.filter_version = 'v_heq_mf10';

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
  '.gs_l_T0.56_L600.v_heq_mf10_m0_d-1_a500.v_heq_mf10';
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
% Algorithm to be used - must match with the function called
config.superpixel_2_seg.method = 'grayscale_ladder';
% Superpixel method: For the 1st stage in superpixel_2_seg this should
% always be config.superpixel_choose.choice.method. So leave the following
% empty ('') for the 1st stage.
config.superpixel_2_seg.superpixel_method = '';
% Superpixel version to be used as basis. For the 1st stage in
% superpixel_2_seg this should always be
% config.superpixel_choose.choice.method. So leave the following empty ('')
% for the 1st stage. 
config.superpixel_2_seg.superpixel_suffix = '';
% thresholds on the boundary field
config.superpixel_2_seg.f_thresholds = 0.58;
% minimum area of a segment - below this they merged with neighbors
config.superpixel_2_seg.area_thresholds = 600;
% Whether mitochondria are to be used to suppress false boundaries
config.superpixel_2_seg.use_mitochondria = false;
% The confidence threshold to be applied on the mitochondria
config.superpixel_2_seg.mitochondria_confidence_threshold = 0.0;
% The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel_2_seg.mitochondria_erosion = 20;
% Whether to use vesicle detection to suppress false boundaries
config.superpixel_2_seg.use_vesicle = false;
% The threshold for the vesicle detection confidence
config.superpixel_2_seg.vesicle_threshold = 0;
% Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.
config.superpixel_2_seg.vesicle_dilation = 0;
% Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel_2_seg.filter_version = 'v_heq_mf10';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9 Choose a segmentation parameter for further processing
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two options:
% --- (a) One parameter for all images ---
% Superpixel method
config.segmentation_choose.choice.method = 'grayscale_ladder';
% use segmentation parameters from previous step
% E.g., for grayscale-ladder it would be of the form
% .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria m<threshold_d<amount of erosion>>  
config.segmentation_choose.choice.seg_suffix = ...
  [ '.gs_l_sp_T0.58_L600', ...
    config.superpixel_choose.choice.seg_suffix, ...
    '.v_heq_mf10'];
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
% (a) 3D annotation file to be used for training
config.linkage.train.manual_annotation_file = ...
  '/groups/chklovskii/medulla/reconstructions/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/2D_segmentation_results/segmentation_stacks_for_linkage_evaluation/seg.wcat.as_2_gt.histeq.region.crop4_global_alignment_0161_0860.296.305.cr928967970826_40.ms3_296.305_3.5k.3.5k.mat';
% '/groups/chklovskii/medulla/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/manual_annotations/region.crop4_global_alignment_0161_0860.296.305.cr928967970826_40/ms3_296.305_3.5k.3.5k/katerina.052609/proofread.katerina.052609.wcat.s.mat';
% () Stack resolution parameters for area and boundary histograms. These
% are to be used if the models having been trained on images of resolution
% different from the current one. This might improve the results,
% especially during boot-strapping.
config.stack.length_factor = 1;
config.stack.area_factor = config.stack.length_factor.^2;
% (c) Type for feature to be used for linkage
config.linkage.feature.type = 'overlap_hist';
config.linkage.feature.bin_size = 10;
% (e) Type of classifier
config.linkage.model.type = 'boost';
% (f) Number of iterations of boosting
config.linkage.model.n_iteration = 40;
% (g) Tree depth
config.linkage.model.tree_depth = 2;
% (b) Version name for the linkage model
config.linkage.train.save_suffix = [ ...
  '.', num2str(config.linkage.model.n_iteration), ...
  '.', num2str(config.linkage.model.tree_depth)];
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
'/groups/chklovskii/medulla/reconstructions/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/2D_segmentation_results/segmentation_stacks_for_linkage_evaluation/seg.v1.s.wcat.region.crop4_global_alignment_0161_0860.296.305.cr928967970826_40.ms3_296.305_3.5k.3.5k.mat';
% '/groups/chklovskii/medulla/reconstructions/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/2D_segmentation_results/segmentation_stacks_for_linkage_evaluation/seg.wcat.region.crop4_global_alignment_0161_0860.296.305.cr928967970826_40.ms3_296.305_3.5k.3.5k.mat';
% Whether to train on the segmentation maps in segmentation_file or on the
% manual_annotation.
config.linkage.evaluate.is_trained_on_manual_annotation = false;
% config.linkage.evaluate.k_fold = 2;
% Linkage thresholds to used for generating the
% false-acceptance/false-rejection curves.
config.linkage.evaluate.link_thresholds = 0:0.1:2;
% whether to display the link decisions
config.linkage.evaluate.display_links = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. Reconstruction Routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pipeline(module_id, config);

return
end
