function config = recon_test_medulla_stack_1_683_686(config, case_ids, is_verbose, is_verbose_figures)
% config = recon_test_medulla_stack_1_683_686(config, case_ids, is_verbose,
%   is_verbose_figures)
%
% Script to test reconstruction software on a small medulla stack.
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  08072009    new protocol for calling modules.
% v1  09102009    test on small medulla stack
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. Constants - Not to be changed by users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Stack root directory redefinition:
% It is recommended that data be stored in
% chklovskiilab/electron_microscopy_data/. The following is for the case
% when the images are stored on an external hard-drive.
config.stack.dir = '/groups/chklovskii/medulla/';
%%% Reconstruction root directory redefinition:
% It is recommended that all intermediate files generated during
% reconstruction be stored in the user's home directory
% em_reconstruction/reconstructions/. However, the directory can be redefined
% to a common share if multiple users wish to read and write on the same
% files and to avoid duplicate storage for large projects.

% Scale for save the coarse al for preliminary viewing and setting ROI.
config.coarse_al.scale = 4; % scale down factor for al_coarse.tif file.

% Whether debug mode is enabled. The program is more verbose.
config.DEBUG = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% II. Configuration parameters - to be set by user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Stack parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
% stack name
config.stack.name = 'test_medulla_stack_1';   
% image name:
% Option 1: If the files are image_0001.tif, image_0002.tif, etc., then
% the image_prefix would be "image_%04d" and the image_suffix would be
% ".tif".
% For files named as a001.tif, a002.tif, etc., the image_prefix would be
% "a%03d" and the image_suffix would be ".tif". 
% If the images are a001.tif, ..., a010.tif then the
% case_ids would be 1:10.
% Option 2: If the images correspond to a trakEM project then
% * specify the trakEM file structure xml file name
config.stack.image_structure = 'test_medulla_0683_0686.xml';
% config.stack.image_structure_differential = 'test_medulla_0683_0684.xml';
% The z-sections within the stack that form the region. If unspecified,
% this defaults to the config.stack.case_ids. Together with
% config.stack.image_structure, config.region.case_ids decides the region
% directory and hence the data_to_be_proofread directory. This should be
% kept constant a particular region's reconstruction but MUST BE CHANGED
% when shifting to a different region.
config.region.case_ids = case_ids;
% * slice ids: specify the z-plane numbers.
if(nargin>1)
  config.stack.case_ids = case_ids;
else
  config.stack.case_ids = 683:684;
end
% ROI: Region Of Interest for the stack. If the entire slice image is
% to be processed then assign empty ([]). Otherwise get_roi(config);
config.stack.roi = [];
% config.stack.roi = get_roi(config);
% Stack resolution parameters for area and boundary histograms. These
% are to be used if the models having been trained on images of resolution
% different from the current one. This might improve the results,
% especially during boot-strapping.
config.stack.length_factor = 1;
config.stack.area_factor = config.stack.length_factor.^2;
% Set the scale at which segmentation should be performed. This is kept
% constant for all superpixel and superpixel-to-seg stages. When stitching
% the segmentations, the label maps can be further downsampled if needed,
% but not artificially upsampled. During linkage and proofreading, the
% segmentations are upsampled to original resolution.
config.stack.segmentation_scale = 2;
config.stack.segmentation_scaling_method = 'bilinear';
% Fold detection and consideration during reconstruction
% Whether folds are considered during segmentation
config.stack.fold.is_considered_in_segmentation = true;
% Whether folds are considered during alignment
config.stack.fold.is_considered_in_alignment = true;

%%%%%%
% 1.1 Precompute data-structures for subsequent routines. E.g., fold masks.
%%%%%%
%%% 1.1.fold
% Minimum area of a component in the thresholded image for it to be
% considered as a fold.
config.precompute.fold.min_fold_area = 1000;
%%% 1.1.BEL
% Whether to compute BEL
config.precompute.boundary.BEL.is_enabled = true; %#### put true to turn on
% Whether to dump test images. Once the test images have been dumped, they
% can be used for multiple parameter versions
config.precompute.boundary.BEL.dump_test_images = true; %true;
% Generate tiles of the images becuase currently (03132009) BEL cannot
% handle 4k x 4k images.
bt(1).minx = 1; bt(1).maxx = 2100; bt(1).miny = 1; bt(1).maxy = 2100;
bt(2).minx = 1; bt(2).maxx = 2100; bt(2).miny = 1950; bt(2).maxy = 4056;
bt(3).minx = 1950; bt(3).maxx = 4056; bt(3).miny = 1950; bt(3).maxy = 4056;
bt(4).minx = 1950; bt(4).maxx = 4056; bt(4).miny = 1; bt(4).maxy = 2100;
config.precompute.boundary.BEL.tiles = bt;
% Whether to print messages
config.precompute.boundary.BEL.is_verbose = is_verbose;
% Version name for the parameters
config.precompute.boundary.BEL.version = 'c_5';
% grey sclae or color images
% 1: color 0: grey scale
config.precompute.boundary.BEL.color = 0;
% Size of the patch around a pixel for representing boundary information
config.precompute.boundary.BEL.patch_size = 35;
% Use pre-computed mask, 0: no mask 1: Canny 2: user defined mask
config.precompute.boundary.BEL.mask = 0;
% use Canny detector at different scales as features
config.precompute.boundary.BEL.Canny=1;
config.precompute.boundary.BEL.set2=0;
% use position information, e.g. row and col, as candidate features, 0: no,
% 1: yes 
config.precompute.boundary.BEL.position=0;
% to stop training on the node below the portion
config.precompute.boundary.BEL.factor_to_stop=0.005;
% maximum depth of PBT
% this largely decides the quality and speed of the final detector
config.precompute.boundary.BEL.depth=30;
% cascade level
config.precompute.boundary.BEL.cascade=0;
% number of weak classifiers in each AdaBoost node
config.precompute.boundary.BEL.number_of_weak_classifiers=50;
% Post processing filter of the BEL results before tiling and saving as
% boundary maps
config.precompute.boundary.BEL.post_process_filter = 'mf7';

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
% Display results using figures?
config.align.precompute.SIFT.is_verbose_figures = is_verbose_figures;
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
% config.align.in_section_align.deformable_mesh.input_tform_format = 'tif_txt';
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
% config.align.in_section_align.SIFT.is_verbose_figures = is_verbose_figures;
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
  '~/temp/sec_align_test_%04d.tif';

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
% config.align.global_align.deformable_mesh.input_tform_format = 'tif_txt';
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
if(~isempty(config.stack.case_ids))
  config.align.global_align.display.save_as_stack_file_name = ...
    [ '~/temp/medulla_align_test_global_', num2str(config.stack.case_ids(1)), ...
    '.', num2str(config.stack.case_ids(end)), '.tif'];
end
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
exceptions = {};
exceptions{1,1} = [685,686]; exceptions{1,2} = '-NCP=200000';
exceptions{2,1} = [686,685]; exceptions{2,2} = '-NCP=200000';
config.align.linkage_align.deformable_mesh.exceptions = exceptions;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% name for the reconstruction
config.reconstruction.name = 'overlap_hist_685.686_083109';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Mitochondria
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
% 5. Vesicle
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.vesicle.apply.dir = '';
config.vesicle.apply.save_suffix = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Superpixel - multi stage
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.superpixel = repmat(config.superpixel, [2 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.1 Superpixel - agglomerative mean boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel(1).method = 'grayscale_AglMeanB';
% () thresholds on the boundary field
config.superpixel(1).f_threshold_seq = 0.01:0.01:0.06;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel(1).save_f_thresholds = 0.06;
% Whether to save the segmentation overlaps as TIF files for viewing.
config.superpixel(1).save_segmentation_overlay = false;
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
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel(1).filter_version = 'BELc\_5\_mf7_e2_d2';
% () whether to print messages
config.superpixel(1).is_verbose = is_verbose;
% () whether to display intermediate results
config.superpixel(1).is_verbose_figures = is_verbose_figures;

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
  s{end+1} = ...
    [ '.gs_amb_T', num2str(f,'%g'), '_', num2str(f-0.01,'%g'), ...
      '_b', num2str(0*(f>0),'%d'), ...
      '.BELc_5_mf7_e2_d2_m0_d5_a40_f', ...
      '.BELc_5_mf7_e2_d2']; %#ok<AGROW>
end
config.superpixel(2).superpixel_suffix = s; 
% () thresholds on the boundary fieldss
config.superpixel(2).f_thresholds = 0.005; %[0.01, 0.05, 0.1, 0.15, 0.2];
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel(2).save_f_thresholds = config.superpixel(2).f_thresholds;
% Whether to save the segmentation overlaps as TIF files for viewing.
config.superpixel(2).save_segmentation_overlay = false;
% () minimum area of a segment - below this they merged with neighbors
config.superpixel(2).area_thresholds = 120;
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel(2).use_mitochondria = false;
% () The confidence threshold to be applied on the mitochondria
config.superpixel(2).mitochondria_confidence_threshold = 0.5; 
% () The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel(2).mitochondria_erosion = 2; % 20
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
config.superpixel(2).filter_version = 'BELc\_5\_mf7_e2_d2';
% () whether to print messages
config.superpixel(2).is_verbose = is_verbose;
% () whether to display intermediate results
config.superpixel(2).is_verbose_figures = is_verbose_figures;

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
    [ '.gs_l_sp_T0.005_L120', ...
      '.gs_amb_T0.06_0.05_b0', ...
      '.BELc_5_mf7_e2_d2_m0_d5_a40_f', ...
      '.BELc_5_mf7_e2_d2', ...
      '.BELc_5_mf7_e2_d2_f'];
% % --- (b) Set a parameter for a batch of sections ---
% % Useful if blocks of sections have different parameters.
% config.superpixel_choose.set_choice.method = 'grayscale_ladder';
% config.superpixel_choose.set_choice.seg_suffix = ...
%   [ '.gs_l_sp_T0.005_L120', ...
%       '.gs_amb_T0.06_0.05_b0', ...
%       '.BELc_5_mf7_e2_d2_m0_d5_a40_f', ...  
%       '.BELc_5_mf7_e2_d2', ...
%       '.BELc_5_mf7_e2_d2_f'];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8 Multistage superpixel to segment computation
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.superpixel_2_seg = repmat(config.superpixel_2_seg, [2 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.1 Superpixel to segment - agglomerative mean boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel_2_seg(1).method = 'grayscale_AglMeanB';
% Superpixel method: For the 1st stage in superpixel_2_seg this should
% always be config.superpixel_choose.choice.method. So leave the following
% empty ('') for the 1st stage.
config.superpixel_2_seg(1).superpixel_method = '';
% (b) Superpixel version to be used as basis. For the 1st stage in
% superpixel_2_seg this should always be
% config.superpixel_choose.choice.method. So leave the following empty ('')
% for the 1st stage. 
config.superpixel_2_seg(1).superpixel_suffix = '';
% (c) thresholds on the boundary field
config.superpixel_2_seg(1).f_threshold_seq = 0.01:0.01:0.11;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel_2_seg(1).save_f_thresholds = 0.11;
% Whether to save the segmentation overlaps as TIF files for viewing.
config.superpixel_2_seg(1).save_segmentation_overlay = false;
% () minimum boundary length (pixel units) for two segments to be merged
% should increase with increasing f_thresholds. E.g., 0 for f_thresholds
% <=0.65 and 20 otherwise
config.superpixel_2_seg(1).length_threshold_seq = ...
  0*(config.superpixel_2_seg(1).f_threshold_seq>0.65);
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel_2_seg(1).use_mitochondria = false;
% () The confidence threshold to be applied on the mitochondria
config.superpixel_2_seg(1).mitochondria_confidence_threshold = 0.0;
% () The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel_2_seg(1).mitochondria_erosion = 20;
% () Whether to use vesicle detection to suppress false boundaries
config.superpixel_2_seg(1).use_vesicle = false;
% () The threshold for the vesicle detection confidence
config.superpixel_2_seg(1).vesicle_threshold = 0;
% () Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.
config.superpixel_2_seg(1).vesicle_dilation = 0;
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel_2_seg(1).filter_version = 'BELc\_5\_mf7_e2_d4';
% () whether to print messages
config.superpixel_2_seg(1).is_verbose = is_verbose;
% () whether to display intermediate results
config.superpixel_2_seg(1).is_verbose_figures = is_verbose_figures;
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.2 Superpixel to segment - agglomerative median boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel_2_seg(2).method = 'grayscale_AglMedianB';
% () Superpixel method
config.superpixel_2_seg(2).superpixel_method = 'grayscale_AglMeanB';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria m<threshold_d<amount of erosion>>
config.superpixel_2_seg(2).superpixel_suffix = ...
  [ '.gs_amb_sp_T0.11_0.1_b0', ...
    config.superpixel_choose.choice.seg_suffix, ...
    '.BELc_5_mf7_e2_d4_f'];
% '.gs_l_T0.47_L100.v_heq_mf3_e2_m0_d2_a500'; 
% (c) thresholds on the boundary field
config.superpixel_2_seg(2).f_threshold_seq = 0.01:0.01:0.08;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel_2_seg(2).save_f_thresholds = 0.08;
% Whether to save the segmentation overlaps as TIF files for viewing.
config.superpixel_2_seg(2).save_segmentation_overlay = false;
% () minimum boundary length (pixel units) for two segments to be merged
% should increase with increasing f_thresholds. E.g., 0 for f_thresholds
% <=0.65 and 20 otherwise
config.superpixel_2_seg(2).length_threshold_seq = ...
  0*(config.superpixel_2_seg(2).f_threshold_seq>0.65);
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel_2_seg(2).use_mitochondria = false;
% () The confidence threshold to be applied on the mitochondria
config.superpixel_2_seg(2).mitochondria_confidence_threshold = 0.0;
% () The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel_2_seg(2).mitochondria_erosion = 20;
% () Whether to use vesicle detection to suppress false boundaries
config.superpixel_2_seg(2).use_vesicle = false;
% () The threshold for the vesicle detection confidence
config.superpixel_2_seg(2).vesicle_threshold = 0;
% () Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.
config.superpixel_2_seg(2).vesicle_dilation = 0;
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel_2_seg(2).filter_version = 'BELc\_5\_mf7_e2_d4';
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
config.segmentation_choose.choice.method = 'grayscale_AglMedianB'; %'grayscale_AglBIFL';
% use segmentation parameters from previous step
% E.g., for grayscale-ladder it would be of the form
% .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria m<threshold_d<amount of erosion>>  
config.segmentation_choose.choice.seg_suffix = ...
  [ '.gs_amdb_sp_T0.08_0.07_b0', ...
    '.gs_amb_sp_T0.11_0.1_b0', ...
    config.superpixel_choose.choice.seg_suffix, ...
    '.BELc_5_mf7_e2_d4_f', ...
    '.BELc_5_mf7_e2_d4_f'];
% --- (b) Set a parameter for a batch of sections ---
% Useful if blocks of sections have different parameters.
config.segmentation_choose.set_choice.method = 'grayscale_AglMedianB';
config.segmentation_choose.set_choice.seg_suffix = ...
  [ '.gs_amdb_sp_T0.08_0.07_b0', ...
    '.gs_amb_sp_T0.11_0.1_b0', ...
    config.superpixel_choose.choice.seg_suffix, ...
    '.BELc_5_mf7_e2_d4_f', ...
    '.BELc_5_mf7_e2_d4_f'];
  
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Align segmentations of tiles within a section and correspond
% overlapping segments
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.align_segmentation.is_verbose = is_verbose;
config.align_segmentation.method = 'boundary_bipartite_match';
config.align_segmentation.delta_margin = 7;
config.align_segmentation.min_n_votes = 15;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11. 3D linkage graph - training
%%%%%%%%%%%%%%%%%%%%%%%%%%
% () Stack resolution parameters for area and boundary histograms. These
% are to be used if the models having been trained on images of resolution
% different from the current one. This might improve the results,
% especially during boot-strapping.
config.stack.length_factor = 1;
config.stack.area_factor = config.stack.length_factor.^2;
% (a) 3D annotation file to be used for training
config.linkage.train.manual_annotation_file = ...
  '/media/FreeAgent_Drive/em_reconstruction/reconstructions/test_medulla_stack_1/linkage_3D/training_data/proofread.shinya.12152008.wcat.mat';
% (b) Version name for the linkage model
config.linkage.train.save_suffix = '.08312009'; % better to change this into the date you do recon
% (c) Type for feature to be used for linkage
config.linkage.feature.type = 'overlap_hist';
config.linkage.feature.bin_size = 10;
% % (c) Type for feature to be used for linkage
% config.linkage.feature.type = 'intensity_pair_hist_v2c';
% % (d) For intensity pair hist., the intensity bins in the histograms
% config.linkage.feature.intensity_bins = 0:0.1:1.0;
% (e) Type of classifier
config.linkage.model.type = 'boost';
% (f) Number of iterations of boosting
config.linkage.model.n_iteration = 15;
% (g) Tree depth
config.linkage.model.tree_depth = 2;
% (h) The model version to be used when applying the linkage on the stack
config.linkage.apply.model_suffix = '.08312009'; % better to change this into the date you do recon
% (j) Verbose
config.linkage.apply.is_verbose = is_verbose;
% Whether to include proofread linkage graphs and if so, the configuration
% to load from. The proofread linkage graphs are stored with suffix
% [ '.proofread.', config.linkage.include_proofread.proofread_name, ...
%   config.linkage.include_proofread.model_suffix];
config.linkage.include_proofread.is_enabled = true;
% proofread name (same as used for config.proofreader.import.proofread_name
config.linkage.include_proofread.proofread_name = ''; %'ms3_296.305_3.5k.3.5k';
% proofread linkage graph version id.
config.linkage.include_proofread.model_suffix = '.v0';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12. Dump to proofreader
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Which proofreader is being used. E.g., matlab_gui, Raveler.
config.proofreader.method = 'Raveler'; %'matlab_gui';
config.proofreader.roi.xmin = 3704; %4160;
config.proofreader.roi.ymin = 3832; %2984;
config.proofreader.roi.xmax = 6703; %7159;
config.proofreader.roi.ymax = 6831; %5983;
config.proofreader.export.is_verbose_figures = is_verbose_figures;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 13. Evaluate reconstruction volume
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.final_volume.link_thresholds = -.1:0.1:2;
config.evaluate_volume.groundtruth_file = 'manual_annotation.1_5.Alex.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 14. Output various datastructures to trakEM
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Whether to print messages
config.trakEM.output_patch_masks_transforms.is_verbose = is_verbose;
% File name of the output javascript with full path
config.trakEM.output_patch_masks_transforms.output_script_name = ...
  '~/temp/test_patch_trakEM.js';

return
end
