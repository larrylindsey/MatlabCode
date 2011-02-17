function reconstructions_Marta_Win_E_20081209_shiv(module_id)
% reconstruction_ex_Alex_volume(module_id)
% Example for using the reconstruction pipeline
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
% v1  12022008  modified for Marta's spanish Windows machine
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
% config.stack.dir = '';
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
config.stack.name = 'lamina.OTO.shinya.leginon.3500x.July2008';   
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
config.stack.image_structure = 'Ecartridge_e.xml';
% * slice ids: specify the z-plane numbers.
config.stack.case_ids = 0:2;
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
config.stack.segmentation_scale = 2;
config.stack.segmentation_scaling_method = 'bilinear';
% () Fold detection and consideration during reconstruction
% Whether folds are considered during segmentation
config.stack.fold.is_considered_in_segmentation = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Alignment parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%
% 2.1 Precompute data-structures for subsequent alignment routines. E.g.,
% for SIFT based alignment, compute SIFT feature points, matches within
% sections and between adjacent sections
%%%%%%
%%% 2.1.SIFT: Feature points
% Configurations for SIFT. Enable/disable SIFT precomputations
config.align.precompute.SIFT.is_enabled = false;
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

% Whether to compute deformable mesh based correspondences between pairs of
% tiles within a section. This may be useful if SIFT fails.
config.align.precompute.deformable_mesh_in_plane_overlap.is_enabled = true;
% Print messages on prompt?
config.align.precompute.deformable_mesh_in_plane_overlap.is_verbose = true;

%%%%%%
% 2.2 For aligning tiles within a section for stitching the segmentations
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
% Plan B: what to do if SIFT fails
% SIFT may be unable to compute correspondences if the area of overlap
% between tiles is too small. We can use deformable mesh-based
% correspondences in place of SIFT.
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

% %%% The following are "reasonable" parameters if using SIFT
% % Method to be used for alignment, e.g., 'SIFT', 'deformable_mesh'
% config.align.linkage_align.method = 'SIFT';
% % Whether to print intermediate messages
% config.align.linkage_align.SIFT.is_verbose = true;
% % Number of iterations of RANSAC. Transformation with minimum error
% % among all iterations is chosen.
% config.align.linkage_align.SIFT.RANSAC.n_iter = 100;
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
config.align.global_align.method = 'SIFT';
% Whether to print intermediate messages
config.align.global_align.SIFT.is_verbose = true;
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
config.align.global_align.SIFT.min_n_matches = 50;
% Consistency threshold for RANSAC to consider aligning a pair of images
config.align.global_align.SIFT.RANSAC.consistency=0.3;
% Whether to show the images inverted for easier assessment of alignment
config.align.global_align.SIFT.is_inverted_display = true;
% Whether to save the aligned stack into a 3D tif.
config.align.global_align.SIFT.save_as_stack_file_name = '';
% Whether to apply a filter on the images before display
config.align.global_align.SIFT.filter_version = '';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. Reconstruction Routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pipeline(module_id, config);

return

end