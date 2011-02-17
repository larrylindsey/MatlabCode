function seg_eval_BEL_medulla_Leginon_681_684_subset(module_id)
% seg_eval_BEL_medulla_Leginon_681_684_subset(module_id)
% Experiment with BEL on medulla data.
%
% module_id
%     0.0,...,0.20      precomputation
%       General stack precomputation
%       0.01    compute fold masks, train BEL.
%       0.02    generate test images for BEL.
%       0.03    copy BEL test results into boundary maps.
%       Alignment precomputation
%       0.01   SIFT feature point detection, deformable mesh, normalized
%               cross correlation.
%       0.02   SIFT within plane matches
%       0.03   SIFT inter-plane matches
%     0.3               global alignment
%     0.31              display global alignment - can also dump a 3D tif
%                         stack
%     0.32              output aligned stack as a MATLAB variable
%     0.6               align tiles within a section
%     0.61              display sectionwise aligned tiles
%     0.8               align pairs tiles in adjacent sections
%     0.81              display tile-pair alignment
%   1     mitochondria manual annotation
%   2     mitochondria collect training samples
%   3     mitochondria train detector
%   4     mitochondria apply detector
%   5     superpixel segmentation. 5.i to run stage ith stage in superpixel
%           segmentation
%   6     superpixel_2_segment_grayscale_ladder
%   7     linkage_3D_train_intensity_pair_boost
%   8     linkage_3D_gen_linkage_gph_intensity_pair_boost
%   9     generate_al
%         generate_cat
%         generate_superpixel_2_seg_map
%         process_superpixel_2_2Dseg_linkage_graph
%   10    prepare_final_reconstruction_volume
%   11    save_reconstruction_info
%         save_reconstruction_param_to_xml
%   12    get_q_score
%   17    in-plane segmentation alignment (stitching)
%   18    choose superpixel parameters using GUI
%   19    choose segmentation parameters using GUI
%   20    launch proofreading GUI
%   21    output datastructures to trakEM
%     21.01     output patch masks and global transforms
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
% v6  03132009  segmentation on BEL boundaries.
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
config.stack.dir = '/groups/chklovskii/medulla/';
%%% Reconstruction root directory redefinition:
% It is recommended that all intermediate files generated during
% reconstruction be stored in the user's home directory
% em_reconstruction/reconstructions/. However, the directory can redefined
% to a common share if multiple users wish to read and write on the same
% files and to avoid duplicate storage for large projects. For instance,
% uncomment the following line to save intermediate reconstruction files in
% the chklovskii/medulla/ directory:
% config.reconstruction.root_dir = '/groups/chklovskii/medulla/reconstructions/';

config.DEBUG = true;


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
config.stack.image_structure = 'test_0681_0684.xml';
% * slice ids: specify the z-plane numbers.
config.stack.case_ids = 681:684;
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
% % () Set the scale at which segmentation should be performed. This is kept
% % constant for all superpixel and superpixel-to-seg stages. When stitching
% % the segmentations, the label maps can be further downsampled if needed,
% % but not artificially upsampled. During linkage and proofreading, the
% % segmentations are upsampled to original resolution.
% config.stack.segmentation_scale = 2;
% config.stack.segmentation_scaling_method = 'bilinear';
% () Fold detection and consideration during reconstruction
% Whether folds are considered during segmentation
config.stack.fold.is_considered_in_segmentation = true;
% Whether folds are considered during alignment
config.stack.fold.is_considered_in_alignment = true;

%%%%%%
% 1.1 Precompute data-structures for subsequent routines. E.g., fold masks.
%%%%%%
%%% 1.1.fold
% Whether to compute and store fold masks for future use
config.precompute.fold.is_enabled = false; %#### put true to turn on
% Whether to save tif files of patches
config.precompute.fold.save_patch_tifs = false;
%%% 1.1.BEL
% Whether to compute BEL
config.precompute.boundary.BEL.is_enabled = true; %#### put true to turn on
% Whether to dump test images. Once the test images have been dumped, they
% can be used for multiple parameter versions
config.precompute.boundary.BEL.dump_test_images = true; % true;
% Generate tiles of the images becuase currently (03132009) BEL cannot
% handle 4k x 4k images.
bt(1).minx = 1; bt(1).maxx = 2100; bt(1).miny = 1; bt(1).maxy = 2100;
bt(2).minx = 1; bt(2).maxx = 2100; bt(2).miny = 1950; bt(2).maxy = 4056;
bt(3).minx = 1950; bt(3).maxx = 4056; bt(3).miny = 1950; bt(3).maxy = 4056;
bt(4).minx = 1950; bt(4).maxx = 4056; bt(4).miny = 1; bt(4).maxy = 2100;
config.precompute.boundary.BEL.tiles = bt;
% Whether to print messages
config.precompute.boundary.BEL.is_verbose = true;
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
% 2.1 Precomputation
% Precompute data-structures for subsequent alignment routines. E.g.,
% for SIFT based alignment, compute SIFT feature points, matches within
% sections and between adjacent sections
%%%%%%%%%
%%%%% 2.1.Normalized Cross Correlation Alignment Within Section
% Align pairs of tiles within section using normalized correlation for
% initil approximation and possibly as plan B.
config.align.precompute.norm_cross_corr_in_plane.is_enabled = false; %####put true to run
% Print messages on prompt?
config.align.precompute.norm_cross_corr_in_plane.is_verbose = true;
% Display results using figures?
config.align.precompute.norm_cross_corr_in_plane.is_verbose_figures = false;
%%%%% 2.1.Normalized Cross Correlation Alignment Across Adjacent Sections
% Align pairs of tiles across adjacent section using normalized correlation
% for initial approximation and possibly as plan B.
config.align.precompute.norm_cross_corr_inter_plane.is_enabled = false; %####put true to run
% Print messages on prompt?
config.align.precompute.norm_cross_corr_inter_plane.is_verbose = true;
% Display results using figures?
config.align.precompute.norm_cross_corr_inter_plane.is_verbose_figures = false;
%%%%% 2.1.SIFT
% Configurations for SIFT. Enable/disable SIFT precomputations
config.align.precompute.SIFT.is_enabled = false; %####put true to run
%%% 2.1.SIFT: Feature points
% Configurations for SIFT. Enable/disable SIFT precomputations
% Print messages on prompt?
config.align.precompute.SIFT.feature_point.is_verbose = true;
% Filter images before feeding to SIFT feature point detector. By default
% no filtering or when defined empty.
config.align.precompute.SIFT.feature_point.filter_version = '';
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
config.align.precompute.deformable_mesh_in_plane_overlap.is_enabled = false; %####put true to run
% Print messages on prompt?
config.align.precompute.deformable_mesh_in_plane_overlap.is_verbose = true;
% Display results using figures?
config.align.precompute.deformable_mesh_in_plane_overlap.is_verbose_figures = false;
%%%%% 2.1.Deformable mesh within plane while considering folds (for global
%%%%% alignment)
% Whether to precompute correspondences using deformable mesh incase SIFT
% fails
config.align.precompute.deformable_mesh_in_plane_overlap_with_fold.is_enabled = false; %####put true to run
% Print messages on prompt?
config.align.precompute.deformable_mesh_in_plane_overlap_with_fold.is_verbose = true;
% Display results using figures?
config.align.precompute.deformable_mesh_in_plane_overlap_with_fold.is_verbose_figures = false;

%%%%%%
% 2.2 For aligning tiles in adjacent sections (joint alignment) for linkage
% for stitching segmentations 
%%%%%%
% Method to be used for alignment, e.g., 'SIFT', 'deformable_mesh'
config.align.in_section_align.method = 'SIFT';
% Whether to print intermediate messages
config.align.in_section_align.SIFT.is_verbose = true;
% Number of iterations of RANSAC. Transformation with minimum error
% among all iterations is chosen.
config.align.in_section_align.SIFT.RANSAC.n_iter = 5;
% Consistency threshold for RANSAC to consider aligning 
% a pair of images
config.align.in_section_align.SIFT.RANSAC.consistency=0.3;
% Whether to show the images inverted for easier assessment of alignment
config.align.in_section_align.SIFT.is_inverted_display = true;
% What to do incase SIFT fails, e.g., due to low overlap
config.align.in_section_align.SIFT.plan_B_method = 'deformable_mesh';

%%%%%%
% 2.3 For aligning a pair of tiles in adjacent sections for linkage
%%%%%%
%%% The following are "reasonable" parameters if using Deformable mesh
% Method to be used for alignment, e.g., 'SIFT', 'deformable_mesh'
config.align.linkage_align.method = 'deformable_mesh';
% Whether to print intermediate messages
config.align.linkage_align.deformable_mesh.is_verbose = true;
% Filter images before feeding to SIFT feature point detector. By default
% no filtering or when defined empty.
config.align.linkage_align.deformable_mesh.filter_version = '';
% Whether to use the overlap obtained from cross-correlation to constrain
% the deformable mesh. This may improve results in certain cases. This may
% be disabled if the deformable mesh is performing coarse alignment
% internally. By default is disabled.
config.align.linkage_align.deformable_mesh.use_precomputed_overlap = true;

% %%% The following are "reasonable" parameters if using SIFT
% % Method to be used for alignment, e.g., 'SIFT', 'deformable_mesh'
% config.align.linkage_align.method = 'SIFT';
% % Whether to print intermediate messages
% config.align.linkage_align.SIFT.is_verbose = true;
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
% 2.4 For aligning all tiles in global coordinates
%%%
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

%%% The following are "reasonable" parameters if using Deformable mesh
% Method to be used for alignment, e.g., 'SIFT', 'deformable_mesh'
config.align.global_align.method = 'deformable_mesh';
% Whether to print intermediate messages
config.align.global_align.deformable_mesh.is_verbose = true;
% Display results using figures?
config.align.global_align.deformable_mesh.is_verbose_figures = false;
% Calculate adjacency matrix from trakEM2 tiling
config.align.global_align.deformable_mesh.trakEM2.getAdjacency=true;
% Optimization method. Options:
% * least_squares solves simple linear equations to minimize the Euclidean
% distance between transformed correspondence points.
% * levenberg_marquadt_rigid imposes soft constraints on the affine parameters
% that favor near rigid transformations.
config.align.global_align.deformable_mesh.optimization.method = 'least_squares';
% config.align.global_align.deformable_mesh.optimization.method = 'levenberg_marquadt_rigid';
% Save intermediate files
config.align.global_align.deformable_mesh.save = false;

%%% For displaying the aligned images
% Whether to show the images inverted for easier assessment of alignment
config.align.global_align.display.is_inverted_display = true;
% Whether to save the aligned stack into a 3D tif.
config.align.global_align.display.save_as_stack_file_name = '';
% Whether to apply a filter on the images before display
config.align.global_align.display.filter_version = '';
% What to do incase SIFT fails, e.g., due to low overlap
config.align.global_align.display.plan_B_method = 'deformable_mesh';


%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% name for the reconstruction
config.reconstruction.name = 'test_pipeline_v3';

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
% 5. Superpixel - multi stage
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.superpixel = repmat(config.superpixel, [2 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.1 Superpixel - agglomerative mean boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel(1).method = 'grayscale_AglMeanB';
% () thresholds on the boundary field
config.superpixel(1).f_threshold_seq = 0.01:0.01:0.09;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel(1).save_f_thresholds = 0.09; %0.1:0.05:0.4;
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
config.superpixel(1).watershed_filter_version = 'BEL_mf7_e6_d4';
% Boundary version
config.superpixel(1).boundary_version = 'c_5_mf7';
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel(1).filter_version = 'BEL_mf7_e6_d4';
% () whether to display intermediate results
config.superpixel(1).is_verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.2 Superpixel - grayscale ladder
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel(2).method = 'grayscale_ladder';
% () Superpixel method
config.superpixel(2).superpixel_method = 'grayscale_AglMeanB';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
s = {};
for f = 0.09 %config.superpixel(1).save_f_thresholds
  s{end+1} = ...
    [ '.gs_amb_T', num2str(f,'%g'), '_', num2str(f-0.01,'%g'), ...
      '_b', num2str(0*(f>0),'%d'), ...
      '.BEL_mf7_e6_d4_m0_d5_a40_f', ...
      '.BEL_mf7_e6_d4']; %#ok<AGROW>
end
config.superpixel(2).superpixel_suffix = s; 
% () thresholds on the boundary field
config.superpixel(2).f_thresholds = 0.005; %[0.01, 0.05, 0.1, 0.15, 0.2];
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel(2).save_f_thresholds = config.superpixel(2).f_thresholds;
% Whether to save the segmentation overlaps as TIF files for viewing.
config.superpixel(2).save_segmentation_overlay = true;
% () minimum area of a segment - below this they merged with neighbors
config.superpixel(2).area_thresholds = 240;
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
config.superpixel(2).filter_version = 'BEL_mf7_e6_d4';
% Boundary version
config.superpixel(2).boundary_version = 'c_5_mf7';
% () whether to display intermediate results
config.superpixel(2).is_verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.3 Superpixel - agglomerative boundary vs. interior values
% Fisher Linear Discriminant
%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % (a) Algorithm to be used - must match with the function called
% % % config.superpixel(3).method = 'grayscale_AglBIFL';
% % % % () Superpixel method
% % % config.superpixel(3).superpixel_method = 'grayscale_ladder';
% % % % (b) Superpixel version to be used as basis
% % % % % use superpixel segmentation parameters from previous step
% % % % For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
% % % s = {};
% % % for f = 0.09; %config.superpixel(1).save_f_thresholds
% % %   s{end+1} = ...
% % %     [ '.gs_l_sp_T0.005_L240', ...
% % %       '.gs_amb_T', num2str(f,'%g'), '_', ...
% % %       num2str(f-0.01,'%g'), '_b', num2str(0*(f>0),'%d'), ...
% % %       '.BEL_mf7_e6_d4_m0_d5_a40', ...
% % %       '.BEL_mf7_e6_d4', ...
% % %       '.BEL_mf7_e6_d4']; %#ok<AGROW>
% % % end
% % % config.superpixel(3).superpixel_suffix = s; 
% % % % (c) thresholds on the boundary field
% % % config.superpixel(3).f_threshold_seq = 0.01;
% % % % () thresholds on the boundary field for which segmentations should be
% % % % saved. Useful for reducing disk space.
% % % config.superpixel(3).save_f_thresholds = 0.01;
% % % % Whether to save the segmentation overlaps as TIF files for viewing.
% % % config.superpixel(3).save_segmentation_overlay = true;
% % % % () minimum boundary length (pixel units) for two segments to be merged
% % % % should increase with increasing f_thresholds. E.g., 0 for f_thresholds
% % % % <=0.65 and 20 otherwise
% % % config.superpixel(3).length_threshold_seq = ...
% % %   20*(config.superpixel(3).f_threshold_seq>0);
% % % % () Maximum area of a segment, above this segments are not considered for
% % % % merging.
% % % config.superpixel(3).max_area_threshold = 2000;
% % % % () Whether mitochondria are to be excluded from image statistics
% % % config.superpixel(3).exclude_mitochondria_stat = true;
% % % % () Whether mitochondria are to be used to suppress false boundaries
% % % config.superpixel(3).use_mitochondria = false;
% % % % () The confidence threshold to be applied on the mitochondria
% % % config.superpixel(3).mitochondria_confidence_threshold = 0.0;
% % % % () The thresholded mitochondria mask is eroded to ensure that correct
% % % % boundaries are not obliterated
% % % config.superpixel(3).mitochondria_erosion = 20;
% % % % * Minimum of a thresholded connected component for it to be considered
% % % % mitochondria
% % % config.superpixel(3).mitochondria_min_area = 40;
% % % % () Whether to use vesicle detection to suppress false boundaries
% % % config.superpixel(3).use_vesicle = false;
% % % % () The threshold for the vesicle detection confidence
% % % config.superpixel(3).vesicle_threshold = 0;
% % % % () Radius of the vesicles. Boundaries within this distance from detected
% % % % vesicle centers are obliterated.
% % % config.superpixel(3).vesicle_dilation = 0;
% % % % () Image filtering before performing segmentation (affects boundary
% % % % detection).
% % % config.superpixel(3).filter_version = 'v0';
% % % % () whether to display intermediate results
% % % config.superpixel(3).is_verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7 Choose a superpixel segmentation parameter for further processing
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two options:
% --- (a) One parameter for all images ---
% Superpixel method
config.superpixel_choose.choice.method = 'grayscale_ladder'; %'grayscale_AglBIFL';
% use superpixel segmentation parameters from previous step
% E.g., for grayscale-ladder it would be of the form
% .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria m<threshold_d<amount of erosion>>  
config.superpixel_choose.choice.seg_suffix = ...
    [ '.gs_l_sp_T0.005_L240', ...
      '.gs_amb_T0.09_0.08_b0', ...
      '.BEL_mf7_e6_d4_m0_d5_a40_f', ...
      '.BEL_mf7_e6_d4', ...
      '.BEL_mf7_e6_d4_f'];
% % % config.superpixel_choose.choice.seg_suffix = ...
% % %     [ '.gs_abif_sp_T0.01_0_b20', ...
% % %       '.gs_l_sp_T0.005_L240', ...
% % %       '.gs_amb_T0.09_0.08_b0', ...
% % %       '.BEL_mf7_e6_d4_m0_d5_a40', ...
% % %       '.BEL_mf7_e6_d4', ...
% % %       '.BEL_mf7_e6_d4', ...
% % %       '.v0'];
% --- (b) For each image choose a parameter from a set ---
% Use a GUI to choose between a set of parameters for each image.
% config.superpixel_choose.choice.method = '#(SP_CHOICE)';
% config.superpixel_choose.choice.seg_suffix = '#(SP_CHOICE)';
% % * This is a list of superpixel parameters. Each element is a struct with
% % fields:
% %   .method       the segmentation method, e.g., 'grayscale_ladder'.
% %   .seg_suffix   parameters for segmentation, similar to those used in
% %                   superpixel_suffix in reconstruction scripts.
% config.superpixel_choose.param(1).method = 'grayscale_AglBIFL';
% config.superpixel_choose.param(1).seg_suffix = ...
%   '.gs_abif_sp_T0.26_0.25_b20.gs_amdb_sp_T0.18_0.17_b0.gs_amb_T0.18_0.17_b0.BEL_e4_d4.BEL_e4_d4.v0';
% 
% config.superpixel_choose.param(2).method = 'grayscale_AglBIFL';
% config.superpixel_choose.param(2).seg_suffix = ...
%   '.gs_abif_sp_T0.26_0.25_b20.gs_amdb_sp_T0.18_0.17_b0.gs_amb_T0.18_0.17_b0.BEL_e4_d4.BEL_e4_d4.v0';
% 
% config.superpixel_choose.param(3).method = 'grayscale_AglBIFL';
% config.superpixel_choose.param(3).seg_suffix = ...
%   '.gs_abif_sp_T0.26_0.25_b20.gs_amdb_sp_T0.18_0.17_b0.gs_amb_T0.18_0.17_b0.BEL_e4_d4.BEL_e4_d4.v0';
% 
% config.superpixel_choose.param(4).method = 'grayscale_AglBIFL';
% config.superpixel_choose.param(4).seg_suffix = ...
%   '.gs_abif_sp_T0.26_0.25_b20.gs_amdb_sp_T0.18_0.17_b0.gs_amb_T0.18_0.17_b0.BEL_e4_d4.BEL_e4_d4.v0';
% % * Desired decision tree: This specifies the manner in which choices are
% % to be presented to the user. Two options:
% % (1) use one of the predefined decision trees (see get_basic_config.m), OR
% % (2) define a custom decision tree: See get_basic_config, Sec. 12 for examples.
% % config.superpixel_choose.decision_tree = config.decision_tree_predef_binary_4_choice;
% config.superpixel_choose.decision_tree = config.decision_tree_predef_binary_4_choice_bias_4;
% % Whether to display the image after histogram equalization
% config.superpixel_choose.is_image_histeq = true;

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
config.superpixel_2_seg(1).f_threshold_seq = 0.01:0.01:0.18;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel_2_seg(1).save_f_thresholds = 0.18; %[0.18, 0.26];
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
% Boundary version
config.superpixel_2_seg(1).boundary_version = 'c_5_mf7';
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel_2_seg(1).filter_version = 'BEL_mf7';
% () whether to display intermediate results
config.superpixel_2_seg(1).is_verbose = true;

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
  [ '.gs_amb_sp_T0.18_0.17_b0', ...
    config.superpixel_choose.choice.seg_suffix, ...
    '.BEL_mf7_f'];
% '.gs_l_T0.47_L100.v_heq_mf3_e2_m0_d2_a500'; 
% (c) thresholds on the boundary field
config.superpixel_2_seg(2).f_threshold_seq = 0.01:0.01:0.18;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel_2_seg(2).save_f_thresholds = 0.18;
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
config.superpixel_2_seg(2).filter_version = 'BEL_mf7';
% Boundary version
config.superpixel_2_seg(2).boundary_version = 'c_5_mf7';
% () whether to display intermediate results
config.superpixel_2_seg(2).is_verbose = true;

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % 8.3 Superpixel to segment - agglomerative boundary vs. interior values
% % % % Fisher Linear Discriminant
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % (a) Algorithm to be used - must match with the function called
% % % config.superpixel_2_seg(3).method = 'grayscale_AglBIFL';
% % % % () Superpixel method
% % % config.superpixel_2_seg(3).superpixel_method = 'grayscale_AglMedianB';
% % % % (b) Superpixel version to be used as basis
% % % % % use superpixel segmentation parameters from previous step
% % % % For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
% % % config.superpixel_2_seg(3).superpixel_suffix = ...
% % %   [ '.gs_amdb_sp_T0.18_0.17_b0', ...
% % %     '.gs_amb_sp_T0.18_0.17_b0', ...
% % %     config.superpixel_choose.choice.seg_suffix, ...
% % %     '.BEL_mf7_f', ...
% % %     '.BEL_mf7_f'];
% % % % (c) thresholds on the boundary field
% % % config.superpixel_2_seg(3).f_threshold_seq = 0.01:0.01:0.2;
% % % % () thresholds on the boundary field for which segmentations should be
% % % % saved. Useful for reducing disk space.
% % % config.superpixel_2_seg(3).save_f_thresholds = 0.26; %0.36:0.04:0.6;
% % % % Whether to save the segmentation overlaps as TIF files for viewing.
% % % config.superpixel_2_seg(3).save_segmentation_overlay = true;
% % % % () minimum boundary length (pixel units) for two segments to be merged
% % % % should increase with increasing f_thresholds. E.g., 0 for f_thresholds
% % % % <=0.65 and 20 otherwise
% % % config.superpixel_2_seg(3).length_threshold_seq = ...
% % %   20*(config.superpixel_2_seg(3).f_threshold_seq>0);
% % % % () Maximum area of a segment, above this segments are not considered for
% % % % merging.
% % % config.superpixel_2_seg(3).max_area_threshold = 2000;
% % % % () Whether mitochondria are to be excluded from image statistics
% % % config.superpixel_2_seg(3).exclude_mitochondria_stat = true;
% % % % () Whether mitochondria are to be used to suppress false boundaries
% % % config.superpixel_2_seg(3).use_mitochondria = false;
% % % % () The confidence threshold to be applied on the mitochondria
% % % config.superpixel_2_seg(3).mitochondria_confidence_threshold = 0.0;
% % % % () The thresholded mitochondria mask is eroded to ensure that correct
% % % % boundaries are not obliterated
% % % config.superpixel_2_seg(3).mitochondria_erosion = 20;
% % % % * Minimum of a thresholded connected component for it to be considered
% % % % mitochondria
% % % config.superpixel_2_seg(3).mitochondria_min_area = 40;
% % % % () Whether to use vesicle detection to suppress false boundaries
% % % config.superpixel_2_seg(3).use_vesicle = false;
% % % % () The threshold for the vesicle detection confidence
% % % config.superpixel_2_seg(3).vesicle_threshold = 0;
% % % % () Radius of the vesicles. Boundaries within this distance from detected
% % % % vesicle centers are obliterated.
% % % config.superpixel_2_seg(3).vesicle_dilation = 0;
% % % % () Image filtering before performing segmentation (affects boundary
% % % % detection).
% % % config.superpixel_2_seg(3).filter_version = 'v0';
% % % % () whether to display intermediate results
% % % config.superpixel_2_seg(3).is_verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9 Choose a segmentation parameter for further processing
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two options:
% --- (a) One parameter for all images ---
% Segmentation method
config.segmentation_choose.choice.method = 'grayscale_AglBIFL';
% use segmentation parameters from previous step
% E.g., for grayscale-ladder it would be of the form
% .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria m<threshold_d<amount of erosion>>  
config.segmentation_choose.choice.seg_suffix = ...
  [ '.gs_abif_sp_T0.26_0.25_b20', ...
    '.gs_amdb_sp_T0.18_0.17_b0', ...
    '.gs_amb_sp_T0.18_0.17_b0', ...
    config.superpixel_choose.choice.seg_suffix, ...
    '.BEL_mf7_f', ...
    '.BEL_mf7_f', ...
    '.v0_f'];
% --- (b) For each image choose a parameter from a set ---
% Use a GUI to choose between a set of parameters for each image.
% % config.segmentation_choose.choice.method = '#(SEG_CHOICE)';
% % config.segmentation_choose.choice.seg_suffix = '#(SEG_CHOICE)';
% % % * This is a list of superpixel parameters. Each element is a struct with
% % % fields:
% % %   .method       the segmentation method, e.g., 'grayscale_ladder'.
% % %   .seg_suffix   parameters for segmentation, similar to those used in
% % %                   superpixel_suffix in reconstruction scripts.
% % config.segmentation_choose.param(1).method = 'grayscale_AglBIFL';
% % config.segmentation_choose.param(1).seg_suffix = ...
% %   '.gs_abif_sp_T0.26_0.25_b20.gs_amdb_sp_T0.18_0.17_b0.gs_amb_T0.18_0.17_b0.BEL_e4_d4.BEL_e4_d4.v0';
% % 
% % config.segmentation_choose.param(2).method = 'grayscale_AglBIFL';
% % config.segmentation_choose.param(2).seg_suffix = ...
% %   '.gs_abif_sp_T0.26_0.25_b20.gs_amdb_sp_T0.26_0.25_b0.gs_amb_T0.26_0.25_b0.BEL_e4_d4.BEL_e4_d4.v0';
% % 
% % config.segmentation_choose.param(3).method = 'grayscale_AglBIFL';
% % config.segmentation_choose.param(3).seg_suffix = ...
% %   '.gs_abif_sp_T0.26_0.25_b20.gs_amdb_sp_T0.18_0.17_b0.gs_amb_T0.18_0.17_b0.BEL_e4_d4.BEL_e4_d4.v0';
% % 
% % config.segmentation_choose.param(4).method = 'grayscale_AglBIFL';
% % config.segmentation_choose.param(4).seg_suffix = ...
% %   '.gs_abif_sp_T0.26_0.25_b20.gs_amdb_sp_T0.18_0.17_b0.gs_amb_T0.18_0.17_b0.BEL_e4_d4.BEL_e4_d4.v0';
% * Desired decision tree: This specifies the manner in which choices are
% to be presented to the user. Two options:
% (1) use one of the predefined decision trees (see get_basic_config.m), OR
% (2) define a custom decision tree: See get_basic_config, Sec. 12 for examples.
config.superpixel_choose.decision_tree = config.decision_tree_predef_binary_4_choice;
% Whether to display the image after histogram equalization
config.segmentation_choose.is_image_histeq = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Align segmentations of tiles within a section and correspond
% overlapping segments
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.align_segmentation.is_verbose = false;
config.align_segmentation.method = 'boundary_bipartite_match';
config.align_segmentation.delta_margin = 7;
config.align_segmentation.min_n_votes = 15;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11. 3D linkage graph - training and generation
%%%%%%%%%%%%%%%%%%%%%%%%%%
% () Stack resolution parameters for area and boundary histograms. These
% are to be used if the models having been trained on images of resolution
% different from the current one. This might improve the results,
% especially during boot-strapping.
config.stack.length_factor = 25/35;
config.stack.area_factor = config.stack.length_factor.^2;
% (a) 3D annotation file to be used for training
config.linkage.train.manual_annotation_file = ...
  '/groups/chklovskii/chklovskiilab/electron_microscopy_data/medulla_HPF_2500x_S3-R1.005_007/manual_annotations/proofread.shinya.07282008.wcat.mat';
% (b) Version name for the linkage model
config.linkage.train.save_suffix = '.09102008';
% (c) Type for feature to be used for linkage
config.linkage.feature.type = 'intensity_pair_hist_v2c';
% (d) For intensity pair hist., the intensity bins in the histograms
config.linkage.feature.intensity_bins = 0:0.1:1.0;
% (e) Type of classifier
config.linkage.model.type = 'boost';
% (f) Number of iterations of boosting
config.linkage.model.n_iteration = 40;
% (g) Tree depth
config.linkage.model.tree_depth = 2;
% (h) The model version to be used when applying the linkage on the stack
config.linkage.apply.model_suffix = '.09102008';
% (j) Verbose
config.linkage.apply.is_verbose = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12. Dump to proofreader
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Which proofreader is being used. E.g., matlab_gui, Raveler.
config.proofreader.method = 'matlab_gui'; %'Raveler';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 13. Evaluate reconstruction volume
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.final_volume.link_thresholds = -.1:0.1:2;
config.evaluate_volume.groundtruth_file = 'manual_annotation.1_5.Alex.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 14. Output various datastructures to trakEM
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Whether to print messages
config.trakEM.output_patch_masks_transforms.is_verbose = true;
% File name of the output javascript with full path
config.trakEM.output_patch_masks_transforms.output_script_name = ...
  '~/temp/test_patch_trakEM.js';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. Reconstruction Routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pipeline(module_id, config);

return

end
