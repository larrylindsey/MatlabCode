function test_SIFT_alignment_JC_030209(module_id, case_ids)

% test_SIFT_alignment_JC_030209(module_id)
% Test SIFT-based alignment on JC (Svoboda Lab.) array tomography data
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
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
% v1  08152008  lamina J recon.
% v2  08252008  test SIFT-based alignment
% v3  10222008  test on array tomography data
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. Constants - Not to be changed by users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config = get_basic_config();

config.stack.dir = '/groups/chklovskii/chklovskiilab/JC_data/';

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
config.stack.name = 'S152_1';   
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
config.stack.image_structure = 'S152_1RGB.xml';
% * slice ids: specify the z-plane numbers.
if(nargin>1)
  config.stack.case_ids = case_ids;
else
  config.stack.case_ids = 26:36;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Alignment parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%
% 2.1 Precomputation
% Precompute data-structures for subsequent alignment routines. E.g.,
% for SIFT based alignment, compute SIFT feature points, matches within
% sections and between adjacent sections
%%%%%%%%%
%%%%% 2.1.Normalized Cross Correlation Alignment Within Section
% Align pairs of tiles within section using normalized correlation for
% initil approximation and possibly as plan B.
config.align.precompute.norm_cross_corr_in_plane.is_enabled = false;
% Print messages on prompt?
config.align.precompute.norm_cross_corr_in_plane.is_verbose = true;
% Display results using figures?
config.align.precompute.norm_cross_corr_in_plane.is_verbose_figures = false;
%%%%% 2.1.Normalized Cross Correlation Alignment Across Adjacent Sections
% Align pairs of tiles across adjacent section using normalized correlation
% for initial approximation and possibly as plan B.
config.align.precompute.norm_cross_corr_inter_plane.is_enabled = false;
% Print messages on prompt?
config.align.precompute.norm_cross_corr_inter_plane.is_verbose = true;
% Display results using figures?
config.align.precompute.norm_cross_corr_inter_plane.is_verbose_figures = false;
%%%%% 2.1.SIFT
% Configurations for SIFT. Enable/disable SIFT precomputations
config.align.precompute.SIFT.is_enabled = true;
%%% 2.1.SIFT: Feature points
% Configurations for SIFT. Enable/disable SIFT precomputations
% Print messages on prompt?
config.align.precompute.SIFT.feature_point.is_verbose = true;
% Filter images before feeding to SIFT feature point detector. By default
% no filtering or when defined empty.
config.align.precompute.SIFT.feature_point.filter_version = 'mul200_mf15';
% Downsample images before computing SIFT feature points
config.align.precompute.SIFT.feature_point.downsampling = 512;
%%% 2.1.SIFT: Within-section SIFT matches
% Print messages on prompt?
config.align.precompute.SIFT.in_section_match.is_verbose = true;
% Calculate adjacency matrix from trakEM2 tiling
config.align.precompute.SIFT.in_section_match.trakEM2.getAdjacency = true;
% Save intermediate files
config.align.precompute.SIFT.in_section_match.save = true;
%%% 2.1.SIFT: Inter-section SIFT matches
% Print messages on prompt?
config.align.precompute.SIFT.inter_section_match.is_verbose = true;
% Calculate adjacency matrix from trakEM2 tiling
config.align.precompute.SIFT.inter_section_match.trakEM2.getAdjacency = true;
% Save intermediate files
config.align.precompute.SIFT.inter_section_match.save = true;
%%%%% 2.1.Deformable mesh within plane while ignoring folds (for stitching)
% Whether to precompute correspondences using deformable mesh incase SIFT
% fails
config.align.precompute.deformable_mesh_in_plane_overlap.is_enabled = false;
% Print messages on prompt?
config.align.precompute.deformable_mesh_in_plane_overlap.is_verbose = true;
% Display results using figures?
config.align.precompute.deformable_mesh_in_plane_overlap.is_verbose_figures = false;
%%%%% 2.1.Deformable mesh within plane while considering folds (for global
%%%%% alignment)
% Whether to precompute correspondences using deformable mesh incase SIFT
% fails
config.align.precompute.deformable_mesh_in_plane_overlap_with_fold.is_enabled = false;
% Print messages on prompt?
config.align.precompute.deformable_mesh_in_plane_overlap_with_fold.is_verbose = true;
% Display results using figures?
config.align.precompute.deformable_mesh_in_plane_overlap_with_fold.is_verbose_figures = false;

%%%%%%
% 2.2 For aligning tiles within sections (joint alignment) for stitching
% segmentation
%%%%%%
% Method to be used for alignment, e.g., 'SIFT', 'deformable_mesh'
config.align.in_section_align.method = 'SIFT';
% Whether to print intermediate messages
config.align.in_section_align.SIFT.is_verbose = true;
% Display results using figures?
config.align.in_section_align.SIFT.is_verbose_figures = false;
% Number of iterations of RANSAC. Transformation with minimum error
% among all iterations is chosen.
config.align.in_section_align.SIFT.RANSAC.n_iter = 5;
% Consistency threshold for RANSAC to consider aligning 
% a pair of images
config.align.in_section_align.SIFT.RANSAC.consistency=0.3;
% allowed discrepancy between consistent matches
config.align.in_section_align.SIFT.RANSAC.max_discrepancy_for_consistency = 5;
% Whether to show the images inverted for easier assessment of alignment
config.align.in_section_align.SIFT.is_inverted_display = true;
% What to do incase SIFT fails, e.g., due to low overlap
config.align.in_section_align.SIFT.plan_B_method = 'deformable_mesh';
% Initial sample correspondences to be used for RANSAC. Two options:
% * Fixed number of initial sample correspondences used for computing an initial
% estimate of model. Should be greater 2 for affine transformations.
% config.align.global_align.SIFT.RANSAC.n_initial_samples = 6;
% * Fraction of inliers to carried over - also decides the fraction of best
% inliers used initially
config.align.in_section_align.SIFT.RANSAC.best_inlier_carryover_frac = 0.15;
% Minimum number of matches between two tiles for their matches to be
% included in the section alignment. For pairs of tiles with marginal
% overlap, the matches are dominated by false correspondences. Therefore,
% the purity of union of matches is improved by eliminating small sets.
config.align.in_section_align.SIFT.min_n_matches = 25;

%%% For displaying the aligned images
% Whether to show the images inverted for easier assessment of alignment
config.align.in_section_align.display.is_inverted_display = false;
% Whether to save the aligned stack into a 3D tif.
config.align.in_section_align.display.save_as_stack_file_name = '';
% Whether to apply a filter on the images before display
config.align.in_section_align.display.filter_version = 'g_mul100';
% Whether to scale the images before display
config.align.in_section_align.display.scale = 4;

%%%%%%
% 2.3 For aligning pairs of tiles within sections.
%%%%%%
% Method to be used for alignment, e.g., 'SIFT', 'deformable_mesh'
config.align.in_section_tile_pair_align.method = 'SIFT';
% Whether to print intermediate messages
config.align.in_section_tile_pair_align.SIFT.is_verbose = true;
% Display results using figures?
config.align.in_section_tile_pair_align.SIFT.is_verbose_figures = false;
% Number of iterations of RANSAC. Transformation with minimum error
% among all iterations is chosen.
config.align.in_section_tile_pair_align.SIFT.RANSAC.n_iter = 5;
% Minimum number of matches between two tiles for their matches to be
% included in the section alignment. For pairs of tiles with marginal
% overlap, the matches are dominated by false correspondences. Therefore,
% the purity of union of matches is improved by eliminating small sets.
config.align.in_section_tile_pair_align.SIFT.min_n_matches = 15;
% matches with deviation under current model within this factor of current
% inliers' mean deviation to be added to inliers
config.align.in_section_tile_pair_align.SIFT.RANSAC.inlier_max_relative_deviation = 0.75;
% Consistency threshold for RANSAC to consider aligning 
% a pair of images
config.align.in_section_tile_pair_align.SIFT.RANSAC.consistency=0.3;
% allowed discrepancy between consistent matches
config.align.in_section_tile_pair_align.SIFT.RANSAC.max_discrepancy_for_consistency = 16;
% Initial sample correspondences to be used for RANSAC. Two options:
% * Fixed number of initial sample correspondences used for computing an initial
% estimate of model. Should be greater 2 for affine transformations.
% config.align.global_align.SIFT.RANSAC.n_initial_samples = 6;
% * Fraction of inliers to carried over - also decides the fraction of best
% inliers used initially
config.align.in_section_tile_pair_align.SIFT.RANSAC.best_inlier_carryover_frac = 0.1;
% Whether to show the images inverted for easier assessment of alignment
config.align.in_section_tile_pair_align.SIFT.is_inverted_display = true;
% What to do incase SIFT fails, e.g., due to low overlap
config.align.in_section_tile_pair_align.SIFT.plan_B_method = '';

%%% For displaying the aligned images
% Whether to show the images inverted for easier assessment of alignment
config.align.in_section_align.display.is_inverted_display = false;
% Whether to save the aligned stack into a 3D tif.
config.align.in_section_align.display.save_as_stack_file_name = '';
% Whether to apply a filter on the images before display
config.align.in_section_align.display.filter_version = 'g_mul100';
% Whether to scale the images before display
config.align.in_section_align.display.scale = 4;

%%%%%%
% 2.3 For aligning a pair of tiles in adjacent sections for linkage
%%%%%%
% %%% The following are "reasonable" parameters if using Deformable mesh
% % Method to be used for alignment, e.g., 'SIFT', 'deformable_mesh'
% config.align.linkage_align.method = 'deformable_mesh';
% % Whether to print intermediate messages
% config.align.linkage_align.deformable_mesh.is_verbose = true;
% % Filter images before feeding to SIFT feature point detector. By default
% % no filtering or when defined empty.
% config.align.linkage_align.deformable_mesh.filter_version = '';

%%% The following are "reasonable" parameters if using SIFT
% Method to be used for alignment, e.g., 'SIFT', 'deformable_mesh'
config.align.linkage_align.method = 'SIFT';
% Whether to print intermediate messages
config.align.linkage_align.SIFT.is_verbose = true;
% Number of iterations of RANSAC. Transformation with minimum error
% among all iterations is chosen.
config.align.linkage_align.SIFT.RANSAC.n_iter = 5;
% Minimum number of matches between two tiles for their matches to be
% included in the section alignment. For pairs of tiles with marginal
% overlap, the matches are dominated by false correspondences. Therefore,
% the purity of union of matches is improved by eliminating small sets.
config.align.linkage_align.SIFT.RANSAC.min_n_matches = 25;
% matches with deviation under current model within this factor of current
% inliers' mean deviation to be added to inliers
config.align.linkage_align.SIFT.RANSAC.inlier_max_relative_deviation = 0.75;
% Consistency threshold for RANSAC to consider aligning 
% a pair of images
config.align.linkage_align.SIFT.RANSAC.consistency = 0.3;
% allowed discrepancy between consistent matches
config.align.linkage_align.SIFT.RANSAC.max_discrepancy_for_consistency = 16;
% Initial sample correspondences to be used for RANSAC. Two options:
% * Fixed number of initial sample correspondences used for computing an initial
% estimate of model. Should be greater 2 for affine transformations.
% config.align.linkage_align.SIFT.RANSAC.n_initial_samples = 4;
% * Fraction of inliers to carried over - also decides the fraction of best
% inliers used initially
config.align.linkage_align.SIFT.RANSAC.best_inlier_carryover_frac = 0.1;
% Whether to show the images inverted for easier assessment of alignment
config.align.linkage_align.SIFT.is_inverted_display = true;

%%%
% 2.4 For aligning all tiles in global coordinates
%%%
%%% The following are "reasonable" parameters if using SIFT
% Method to be used for alignment, e.g., 'SIFT', 'deformable_mesh'
config.align.global_align.method = 'SIFT';
% Whether to print intermediate messages
config.align.global_align.SIFT.is_verbose = true;
% Display results using figures?
config.align.global_align.SIFT.is_verbose_figures = false;
% Calculate adjacency matrix from trakEM2 tiling
config.align.global_align.SIFT.trakEM2.getAdjacency=true;
% Save intermediate files
config.align.global_align.SIFT.save = true;
% Number of iterations of RANSAC. Transformation with minimum error
% among all iterations is chosen.
config.align.global_align.SIFT.RANSAC.n_iter = 5;
% Initial sample correspondences to be used for RANSAC. Two options:
% * Fixed number of initial sample correspondences used for computing an initial
% estimate of model. Should be greater 2 for affine transformations.
% config.align.global_align.SIFT.RANSAC.n_initial_samples = 6;
% * Fraction of inliers to carried over - also decides the fraction of best
% inliers used initially
config.align.global_align.SIFT.RANSAC.best_inlier_carryover_frac = 0.25;
% Minimum number of matches between two tiles for their matches to be
% included in the section alignment. For pairs of tiles with marginal
% overlap, the matches are dominated by false correspondences. Therefore,
% the purity of union of matches is improved by eliminating small sets.
config.align.global_align.SIFT.min_n_matches = 25;
% Consistency threshold for RANSAC to consider aligning a pair of images
config.align.global_align.SIFT.RANSAC.consistency=0.3;
% What to do incase SIFT fails, e.g., due to low overlap
config.align.global_align.display.plan_B_method = '';

% %%% The following are "reasonable" parameters if using Deformable mesh
% % Method to be used for alignment, e.g., 'SIFT', 'deformable_mesh'
% config.align.global_align.method = 'deformable_mesh';
% % Whether to print intermediate messages
% config.align.global_align.deformable_mesh.is_verbose = true;
% % Display results using figures?
% config.align.global_align.deformable_mesh.is_verbose_figures = false;
% % Calculate adjacency matrix from trakEM2 tiling
% config.align.global_align.deformable_mesh.trakEM2.getAdjacency=true;
% % Optimization method. Options:
% % * least_squares solves simple linear equations to minimize the Euclidean
% % distance between transformed correspondence points.
% % * levenberg_marquadt_rigid imposes soft constraints on the affine parameters
% % that favor near rigid transformations.
% % config.align.global_align.deformable_mesh.optimization.method = 'least_squares';
% config.align.global_align.deformable_mesh.optimization.method = 'levenberg_marquadt_rigid';
% % Save intermediate files
% config.align.global_align.deformable_mesh.save = true;

%%% For displaying the aligned images
% Whether to show the images inverted for easier assessment of alignment
config.align.global_align.display.is_inverted_display = false;
% Whether to save the aligned stack into a 3D tif.
config.align.global_align.display.save_as_stack_file_name = ...
  ['~/temp/JC_S152_1_', num2str(config.stack.case_ids(1)), '_', ...
  num2str(config.stack.case_ids(end)), '.tif'];
% Whether to apply a filter on the images before display
config.align.global_align.display.filter_version = 'g_mul100';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. Reconstruction Routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pipeline(module_id, config);

return

end
