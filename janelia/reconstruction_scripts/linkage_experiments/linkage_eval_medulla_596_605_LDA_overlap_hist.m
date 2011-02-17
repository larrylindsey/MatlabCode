function config = linkage_eval_medulla_596_605_LDA_overlap_hist(config, case_ids, is_verbose, is_verbose_figures)
% seg_eval_medulla_Leginon_misc_LDA_shrink_LS(module_id)
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
% Option 2: If the images correspond to a trakEM project then
% * specify the trakEM file structure xml file name
config.stack.image_structure = 'eval_data_596_605_1tile.xml';
% config.stack.image_structure = 'crop4_global_alignment_0161_0860.xml';
% The z-sections within the stack that form the region. If unspecified,
% this defaults to the config.stack.case_ids. Together with
% config.stack.image_structure, config.region.case_ids decides the region
% directory and hence the data_to_be_proofread directory. This should be
% kept constant a particular region's reconstruction but MUST BE CHANGED
% when shifting to a different region.
config.region.case_ids = 596:605;
% * slice ids: specify the z-plane numbers.
if(nargin>1)
  config.stack.case_ids = case_ids;
else
  config.stack.case_ids = 596:605;
end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) name for the reconstruction
config.reconstruction.name = 'link_eval';

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
% 2. Alignment parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%
% 2.0 Precomputation
% Precompute data-structures for subsequent alignment routines. E.g.,
% for SIFT based alignment, compute SIFT feature points, matches within
% sections and between adjacent sections
%%%%%%%%%
%%%%% 2.0.SIFT
% Configurations for SIFT. Enable/disable SIFT precomputations
% Print messages on prompt?
config.align.precompute.SIFT.feature_point.is_verbose = is_verbose;
% Filter images before feeding to SIFT feature point detector. By default
% no filtering or when defined empty.
config.align.precompute.SIFT.feature_point.filter_version = '';
% Downsample images before computing SIFT feature points
config.align.precompute.SIFT.feature_point.downsampling = 512;

%%%%% 2.0.Normalized Cross Correlation Alignment Within Section
% Align pairs of tiles within section using normalized correlation for
% initil approximation and possibly as plan B.
% Print messages on prompt?
config.align.precompute.norm_cross_corr_in_plane.is_verbose = is_verbose;
% Display results using figures?
config.align.precompute.norm_cross_corr_in_plane.is_verbose_figures = is_verbose_figures;
%%%%% 2.0.Normalized Cross Correlation Alignment Across Adjacent Sections
% Align pairs of tiles across adjacent section using normalized correlation
% for initial approximation and possibly as plan B.
% Print messages on prompt?
config.align.precompute.norm_cross_corr_inter_plane.is_verbose = is_verbose;
% Display results using figures?
config.align.precompute.norm_cross_corr_inter_plane.is_verbose_figures = is_verbose_figures;

%%%%% 2.0.Deformable mesh within plane while ignoring folds (for stitching)
% Print messages on prompt?
config.align.precompute.deformable_mesh_in_plane_overlap.is_verbose = is_verbose;
% Display results using figures?
config.align.precompute.deformable_mesh_in_plane_overlap.is_verbose_figures = is_verbose_figures;

%%%%%%
% 2.1 For aligning pairs of tiles within a section.
%%%%%%              
%%% The following are "reasonable" parameters if using Deformable mesh
% Method to be used, e.g., 'SIFT', 'deformable_mesh'
config.align.in_section_tile_pair_align.method = 'deformable_mesh';
% Print messages on prompt?
config.align.in_section_tile_pair_align.deformable_mesh.is_verbose = is_verbose;
% Display results using figures?
config.align.in_section_tile_pair_align.deformable_mesh.is_verbose_figures = is_verbose_figures;
% Parameters for deformable mesh
config.align.in_section_tile_pair_align.deformable_mesh.params = '-DFT=0.305';

% %%% The following are "reasonable" parameters if using SIFT
% % Method to be used, e.g., 'SIFT', 'deformable_mesh'
% config.align.in_section_tile_pair_align.method = 'SIFT';
% % Print messages on prompt?
% config.align.in_section_tile_pair_align.SIFT.is_verbose = is_verbose;
% % Calculate adjacency matrix from trakEM2 tiling
% config.align.in_section_tile_pair_align.SIFT.trakEM2.getAdjacency = true;
% % Save intermediate files
% config.align.in_section_tile_pair_align.SIFT.save = true;
% % Number of iterations of RANSAC. Transformation with minimum error
% % among all iterations is chosen.
% config.align.in_section_tile_pair_align.SIFT.RANSAC.n_iter = 5;
% % Consistency threshold for RANSAC to consider aligning 
% % a pair of images
% config.align.in_section_tile_pair_align.SIFT.RANSAC.consistency=0.3;
% % Whether to show the images inverted for easier assessment of alignment
% config.align.in_section_tile_pair_align.SIFT.is_inverted_display = true;
% % What to do incase SIFT fails, e.g., due to low overlap
% config.align.in_section_tile_pair_align.SIFT.plan_B_method = 'deformable_mesh';

%%%%%%
% 2.2 For aligning a pair of tiles in adjacent sections for linkage
%%%%%%
%%% The following are "reasonable" parameters if using Deformable mesh
% Method to be used for alignment, e.g., 'SIFT', 'deformable_mesh'
config.align.linkage_align.method = 'deformable_mesh';
% Whether to print intermediate messages
config.align.linkage_align.deformable_mesh.is_verbose = is_verbose;
config.align.linkage_align.deformable_mesh.is_verbose_figures = is_verbose_figures;
% Filter images before feeding to SIFT feature point detector. By default
% no filtering or when defined empty.
config.align.linkage_align.deformable_mesh.filter_version = '';
% For deformable mesh, whether to read files from matlab_mat or tif_txt (as
% generated by binaries).
config.align.linkage_align.deformable_mesh.input_tform_format = 'tif_txt';
% Parameters for deformable mesh
config.align.linkage_align.deformable_mesh.params = '-DFT=0.305';
% Whether to use the overlap obtained from cross-correlation to constrain
% the deformable mesh. This may improve results in certain cases. This may
% be disabled if the deformable mesh is performing coarse alignment
% internally. By default is disabled.
config.align.linkage_align.deformable_mesh.use_precomputed_overlap = false;

% %%% The following are "reasonable" parameters if using SIFT
% % Method to be used for alignment, e.g., 'SIFT', 'deformable_mesh'
% config.align.linkage_align.method = 'SIFT';
% % Calculate adjacency matrix from trakEM2 tiling
% config.align.linkage_align.SIFT.trakEM2.getAdjacency = true;
% % Save intermediate files
% config.align.linkage_align.SIFT.save = true;
% % Whether to print intermediate messages
% config.align.linkage_align.SIFT.is_verbose = is_verbose;
% % Number of iterations of RANSAC. Transformation with minimum error
% % among all iterations is chosen.
% config.align.linkage_align.SIFT.RANSAC.n_iter = 5;
% % Initial sample correspondences to be used for RANSAC. Two options:
% % * Fixed number of initial sample correspondences used for computing an initial
% % estimate of model. Should be greater 2 for affine transformations.
% % config.align.linkage_align.SIFT.RANSAC.n_initial_samples = 6;
% % * Fraction of inliers to carried over - also decides the fraction of best
% % inliers used initially
% config.align.linkage_align.SIFT.RANSAC.best_inlier_carryover_frac = 0.15;
% % Minimum number of matches between two tiles for their matches to be
% % included in the section alignment. For pairs of tiles with marginal
% % overlap, the matches are dominated by false correspondences. Therefore,
% % the purity of union of matches is improved by eliminating small sets.
% config.align.linkage_align.SIFT.RANSAC.min_n_matches = 50;
% % Consistency threshold for RANSAC to consider aligning 
% % a pair of images
% config.align.linkage_align.SIFT.RANSAC.consistency = 0.3;
% % Whether to show the images inverted for easier assessment of alignment
% config.align.linkage_align.SIFT.is_inverted_display = true;

%%%
% 2.3 For aligning all tiles within a section
%%%
%%% The following are "reasonable" parameters if using Deformable mesh
% Method to be used for alignment, e.g., 'SIFT', 'deformable_mesh'
config.align.in_section_align.method = 'deformable_mesh';
% Whether to print intermediate messages
config.align.in_section_align.deformable_mesh.is_verbose = is_verbose;
% Display results using figures?
config.align.in_section_align.deformable_mesh.is_verbose_figures = is_verbose_figures;
% Calculate adjacency matrix from trakEM2 tiling
config.align.in_section_align.deformable_mesh.trakEM2.getAdjacency=true;

% Optimization method. Options:
% * least_squares solves simple linear equations to minimize the Euclidean
% distance between transformed correspondence points.
% * levenberg_marquadt_rigid imposes soft constraints on the affine parameters
% that favor near rigid transformations.
config.align.in_section_align.deformable_mesh.optimization.method = 'least_squares';
% config.align.global_align.deformable_mesh.optimization.method = 'levenberg_marquadt_rigid';
% For deformable mesh, whether to read files from matlab_mat or tif_txt (as
% generated by binaries).
config.align.in_section_align.deformable_mesh.input_tform_format = 'tif_txt';
% Save intermediate files
config.align.in_section_align.deformable_mesh.save = true;
% What to do incase deformable mesh fails, e.g., due to low overlap.
% 'deformable_mesh': use deformable_mesh_in_plane_overlap - this has
% rectangular patches and is specifically designed for narrow long strips
% of overlap.
config.align.in_section_align.deformable_mesh.plan_B_method = 'deformable_mesh';

% %%% The following are "reasonable" parameters if using SIFT
% % Method to be used for alignment, e.g., 'SIFT', 'deformable_mesh'
% config.align.in_section_align.method = 'SIFT';
% % Whether to print intermediate messages
% config.align.in_section_align.SIFT.is_verbose = true;
% % Display results using figures?
% config.align.in_section_align.SIFT.is_verbose_figures = true;
% % Number of iterations of RANSAC. Transformation with minimum error
% % among all iterations is chosen.
% config.align.in_section_align.SIFT.RANSAC.n_iter = 10;
% % Consistency threshold for RANSAC to consider aligning 
% % a pair of images
% config.align.in_section_align.SIFT.RANSAC.consistency=0.3;
% % allowed discrepancy between consistent matches
% config.align.in_section_align.SIFT.RANSAC.max_discrepancy_for_consistency = 5;
% % Whether to show the images inverted for easier assessment of alignment
% config.align.in_section_align.SIFT.is_inverted_display = true;
% % What to do incase SIFT fails, e.g., due to low overlap
% config.align.in_section_align.SIFT.plan_B_method = 'deformable_mesh';

%%% For displaying the aligned images
% Whether to show the images inverted for easier assessment of alignment
config.align.in_section_align.display.is_inverted_display = true;
% Whether to apply a filter on the images before display
config.align.in_section_align.display.filter_version = '';
% Scale images before display
config.align.in_section_align.display.display_scale = 8;
% Whether to save the aligned sections into a tif.
config.align.in_section_align.display.save_as_tif_prefix = ...
  '/groups/visitors/home/takemuras/temp/sec_align_test_%04d.tif';

%%%
% 2.4 For aligning all tiles in global coordinates
%%%
%%% The following are "reasonable" parameters if using Deformable mesh
% Method to be used for alignment, e.g., 'SIFT', 'deformable_mesh'
config.align.global_align.method = 'deformable_mesh';
% Whether to print intermediate messages
config.align.global_align.deformable_mesh.is_verbose = is_verbose;
% Display results using figures?
config.align.global_align.deformable_mesh.is_verbose_figures = is_verbose_figures;
% Calculate adjacency matrix from trakEM2 tiling
config.align.global_align.deformable_mesh.trakEM2.getAdjacency=true;
% For deformable mesh, whether to read files from matlab_mat or tif_txt (as
% generated by binaries).
config.align.global_align.deformable_mesh.input_tform_format = 'tif_txt';
% Optimization method. Options:
% * least_squares solves simple linear equations to minimize the Euclidean
% distance between transformed correspondence points.
% * levenberg_marquadt_rigid imposes soft constraints on the affine parameters
% that favor near rigid transformations.
 config.align.global_align.deformable_mesh.optimization.method = 'least_squares_rigidified';
% Save intermediate files
config.align.global_align.deformable_mesh.save = true;
% What to do incase deformable mesh fails, e.g., due to low overlap.
% 'deformable_mesh': use deformable_mesh_in_plane_overlap - this has
% rectangular patches and is specifically designed for narrow long strips
% of overlap.
config.align.global_align.deformable_mesh.plan_B_method = 'deformable_mesh';
% Minimum number of correspondences
config.align.global_align.deformable_mesh.min_n_matches = 5;
% Format to save/load global patchwise affine transforms [default '.mat']
config.align.global_align.storage_format = 'txt';

% %Exceptions
% exceptions = {};
% exceptions{1,1} = 116:124; exceptions{1,2} = '-plan_B_method=[]';
% config.align.global_align.deformable_mesh.exceptions = exceptions;

% %%% The following are "reasonable" parameters if using SIFT
% % Method to be used for alignment, e.g., 'SIFT', 'deformable_mesh'
% config.align.global_align.method = 'SIFT';
% % Whether to print intermediate messages
% config.align.global_align.SIFT.is_verbose = true;
% % Calculate adjacency matrix from trakEM2 tiling
% config.align.global_align.SIFT.trakEM2.getAdjacency=true;
% % Save intermediate files
% config.align.global_align.SIFT.save = true;
% % Number of iterations of RANSAC. Transformation with minimum error
% % among all iterations is chosen.
% config.align.global_align.SIFT.RANSAC.n_iter = 5;
% % Initial sample correspondences to be used for RANSAC. Two options:
% % * Fixed number of initial sample correspondences used for computing an initial
% % estimate of model. Should be greater 2 for affine transformations.
% % config.align.global_align.SIFT.RANSAC.n_initial_samples = 6;
% % * Fraction of inliers to carried over - also decides the fraction of best
% % inliers used initially
% config.align.global_align.SIFT.RANSAC.best_inlier_carryover_frac = 0.25;
% % Minimum number of matches between two tiles for their matches to be
% % included in the section alignment. For pairs of tiles with marginal
% % overlap, the matches are dominated by false correspondences. Therefore,
% % the purity of union of matches is improved by eliminating small sets.
% config.align.global_align.SIFT.min_n_matches = 25;
% % Consistency threshold for RANSAC to consider aligning a pair of images
% config.align.global_align.SIFT.RANSAC.consistency=0.3;
% % What to do incase SIFT fails, e.g., due to low overlap
% config.align.global_align.SIFT.plan_B_method = 'deformable_mesh';

%%% For displaying the aligned images
% Whether to show the images inverted for easier assessment of alignment
config.align.global_align.display.is_inverted_display = true;
% Whether to save the aligned stack into a 3D tif.
config.align.global_align.display.save_as_stack_file_name = ...
  [ '/groups/visitors/home/takemuras/temp/medulla_align_test_global_', num2str(config.stack.case_ids(1)), ...
    '.', num2str(config.stack.case_ids(end)), '.tif'];
% Whether to apply a filter on the images before display
config.align.global_align.display.filter_version = '';
% Scale images before display
config.align.global_align.display.display_scale = 8;

% %Exceptions
% %parameters for in-section dmesh align
% config.align.dmesh.intra_section.params = '-DFT=0.305';
% % parameters for across-section dmesh align
% config.align.dmesh.inter_section.params = '-DFT=0.305';
% 
% %Allow for exceptions to be passed to dmesh code
% exceptions = {};
% exceptions{1,1} = [76,77]; exceptions{1,2} = '-SCALE=0.90';
% config.align.linkage_align.deformable_mesh.exceptions=exceptions;

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
config.superpixel(1).watershed_filter_version = 'db_nwo700_cbf2_LDA7_mf7_ps1.5\_0.2_abf_neg';
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel(1).filter_version = 'db_nwo700_cbf2_LDA7_mf7_ps1.5\_0.2_abf_neg';
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
      '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg_m0_d5_a40', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg']; %#ok<AGROW>
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
config.superpixel(2).filter_version = 'db_nwo700_cbf2_LDA7_mf7_ps1.5\_0.2_abf_neg';
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
      '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg_m0_d5_a40', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg']; %#ok<AGROW>
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
config.superpixel(3).filter_version = 'db_nwo700_cbf2_LDA7_mf7_ps1.5\_0.2_abf_neg';
% () whether to print messages
config.superpixel(3).is_verbose = is_verbose;
% () whether to display intermediate results
config.superpixel(3).is_verbose_figures = is_verbose_figures;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.1 Superpixel - prune-classify
%%%%%%%%%%%%%%%%%%%%%%%%%%
if(length(config.superpixel_2_seg)==1)
  config.superpixel_2_seg = repmat(config.superpixel_2_seg, [2 1]);
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
    [ '.gs_l_sp_T1e-05_L300', ...
      '.gs_amdb_sp_T', num2str(f,'%g'), '_', num2str(f-0.01,'%g'), ...
      '_b', num2str(5*(f>-1),'%d'), ...
      '.gs_amb_T0.05_0.04_b0', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg_m0_d5_a40', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg']; %#ok<AGROW>
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
config.superpixel_2_seg(1).filter_version = 'db_nwo700_cbf2_LDA7_mf7_ps1.5\_0.2_abf_neg';
% () whether to print messages
config.superpixel_2_seg(1).is_verbose = is_verbose;
% () whether to display intermediate results
config.superpixel_2_seg(1).is_verbose_figures = is_verbose_figures;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.2 Superpixel - prune-classify
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel_2_seg(2).method = 'shrink_levelset';
% () Superpixel method
config.superpixel_2_seg(2).superpixel_method = 'prune_classify';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
s = {};
for f = config.superpixel(2).save_f_thresholds(4) %0.1
  s{end+1} = ...
    [ '.apc_sp_region_max_2_0.2', ...
      '.gs_l_sp_T1e-05_L300', ...
      '.gs_amdb_sp_T', num2str(f,'%g'), '_', num2str(f-0.01,'%g'), ...
      '_b', num2str(5*(f>-1),'%d'), ...
      '.gs_amb_T0.05_0.04_b0', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg_m0_d5_a40', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
      '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg']; %#ok<AGROW>
end
config.superpixel_2_seg(2).superpixel_suffix = s; 
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
config.superpixel_2_seg(2).filter_version = 'db_nwo700_cbf2_LDA7_mf7_ps1.5\_0.2_abf_neg';
% () whether to print messages
config.superpixel_2_seg(2).is_verbose = is_verbose;
% () whether to display intermediate results
config.superpixel_2_seg(2).is_verbose_figures = is_verbose_figures;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9 Choose a segmentation parameter for further processing
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two options:
% --- (a) One parameter for all images ---
% Segmentation method
% config.segmentation_choose.choice.method = 'grayscale_AglMedianB'; %'grayscale_AglBIFL';
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
  [ '.apc_sp_region_max_2_0.2', ...
    '.gs_l_sp_T1e-05_L300', ...
    '.gs_amdb_sp_T0.4_0.39_b5', ...
    '.gs_amb_T0.05_0.04_b0', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg_m0_d5_a40', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg'];
% --- (b) For each image choose a parameter from a set ---
% % Use a GUI to choose between a set of parameters for each image.
% % * This is a list of superpixel parameters. Each element is a struct with
% % fields:
% %   .method       the segmentation method, e.g., 'grayscale_ladder'.
% %   .seg_suffix   parameters for segmentation, similar to those used in
% %                   superpixel_suffix in reconstruction scripts.
% config.superpixel_choose.param(1).method = 'grayscale_AglBIFL';
% config.superpixel_choose.param(1).seg_suffix = ...
%   '.gs_abif_sp_T0.18_0.17_b2.gs_l_sp_T0.2_L60.gs_amb_T0.63_0.62_b0.v_heq_o1_m0_d5_a40_f.v_heq2.v_heq2_f.v0_f';
% 
% config.superpixel_choose.param(2).method = 'grayscale_AglBIFL';
% config.superpixel_choose.param(2).seg_suffix = ...
%   '.gs_abif_sp_T0.16_0.15_b2.gs_l_sp_T0.2_L60.gs_amb_T0.63_0.62_b0.v_heq_o1_m0_d5_a40_f.v_heq2.v_heq2_f.v0_f';
% 
% config.superpixel_choose.param(3).method = 'grayscale_AglBIFL';
% config.superpixel_choose.param(3).seg_suffix = ...
%   '.gs_abif_sp_T0.14_0.13_b2.gs_l_sp_T0.2_L60.gs_amb_T0.63_0.62_b0.v_heq_o1_m0_d5_a40_f.v_heq2.v_heq2_f.v0_f';
% 
% config.superpixel_choose.param(4).method = 'grayscale_AglBIFL';
% config.superpixel_choose.param(4).seg_suffix = ...
%   '.gs_abif_sp_T0.12_0.11_b2.gs_l_sp_T0.2_L60.gs_amb_T0.63_0.62_b0.v_heq_o1_m0_d5_a40_f.v_heq2.v_heq2_f.v0_f';
% % * Desired decision tree: This specifies the manner in which choices are
% % to be presented to the user. Two options:
% % (1) use one of the predefined decision trees (see get_basic_config.m), OR
% % (2) define a custom decision tree: See get_basic_config, Sec. 12 for examples.
% config.superpixel_choose.decision_tree = config.decision_tree_predef_binary_4_choice;
% % Whether to display the image after histogram equalization
% config.superpixel_choose.is_image_histeq = true;
% --- (c) Set a parameter for a batch of sections ---
% Useful if blocks of sections have different parameters.
% config.segmentation_choose.set_choice.method = 'grayscale_AglMedianB';
% config.segmentation_choose.set_choice.seg_suffix = ...
%   [ '.gs_amdb_sp_T0.02_0.01_b0', ...
%     '.gs_amb_sp_T0.02_0.01_b0', ...
%     config.superpixel_choose.choice.seg_suffix, ...
%     '.BELc_5_mf7_e2_d4_f', ...
%     '.BELc_5_mf7_e2_d4_f'];

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
config.eval_segment_boundary.method = 'prune_classify'; % 'grayscale_ladder';
config.eval_segment_boundary.seg_suffix_format = ...
[ '.apc_sp_region_max_2_0.2', ...
  '.gs_l_sp_T1e-05_L300', ...
  '.gs_amdb_sp%s', ...
  '.gs_amb_T0.05_0.04_b0', ...
  '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg_m0_d5_a40', ...
  '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
  '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
  '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
  '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg'];
% [ '.gs_l_sp_T1e-05_L300', ...
%   '.gs_amdb_sp%s', ...
%   '.gs_amb_T0.05_0.04_b0', ...
%   '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg_m0_d5_a40', ...
%   '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
%   '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
%   '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg'];
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
config.linkage.evaluate.segment_method_automatic = 'prune_classify';
config.linkage.evaluate.segment_suffix_automatic = ...
  [ '.apc_sp_region_max_2_0.2', ...
    '.gs_l_sp_T1e-05_L300', ...
    '.gs_amdb_sp_T0.4_0.39_b5', ...
    '.gs_amb_T0.05_0.04_b0', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg_m0_d5_a40', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1.5_0.2_abf_neg'];
% Automatic linkage suffix
config.linkage.evaluate.linkage_suffix_automatic = ...
  ['.boost.overlap_hist.40.2', config.linkage.evaluate.segment_suffix_automatic];
% Ground-truth segmentation maps to be used for testing.
config.linkage.evaluate.superpixel_method_ground_truth = 'proofread_result';
config.linkage.evaluate.superpixel_suffix_ground_truth = '.prfrd_asp.satoko.042709';
config.linkage.evaluate.segment_method_ground_truth = 'proofread_result';
config.linkage.evaluate.segment_suffix_ground_truth = '.prfrd_seg.satoko.042709';
% Ground-truth linkage suffix
config.linkage.evaluate.linkage_suffix_ground_truth = ...
  '.proofread.satoko.042709.v0#(SEG_CHOICE)';
% Linkage thresholds to used for generating the
% false-acceptance/false-rejection curves.
config.linkage.evaluate.linkage_thresholds = 0.5; % 0:0.1:2;

return
end
