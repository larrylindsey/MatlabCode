function config = recon_medulla_596_605_LDA_seeded_watershed_nomito(config, case_ids, is_verbose, is_verbose_figures)
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
% config.stack.image_prefix = 'a.%03d';
% config.stack.image_suffix = '.tif';
% Option 2: If the images correspond to a trakEM project then
% * specify the trakEM file structure xml file name
% config.stack.image_structure = 'eval_data_596_605_1tile.xml';
config.stack.image_structure = 'crop4_global_alignment_0161_0860.xml';
% The z-sections within the stack that form the region. If unspecified,
% this defaults to the config.stack.case_ids. Together with
% config.stack.image_structure, config.region.case_ids decides the region
% directory and hence the data_to_be_proofread directory. This should be
% kept constant a particular region's reconstruction but MUST BE CHANGED
% when shifting to a different region.
config.region.case_ids = 596:605;
if(nargin>1 && ~isempty(case_ids))
  config.stack.case_ids = case_ids;
else
  config.stack.case_ids = 596:605;
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
config.stack.fold.is_considered_in_segmentation = true;
% Whether folds are considered during alignment
config.stack.fold.is_considered_in_alignment = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) name for the reconstruction
config.reconstruction.name = 'eval_sws_nm_lp_2k';

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
% Whether to show only specific tiles
config.align.linkage_align.display.tile_1 = 5;
config.align.linkage_align.display.tile_2 = 9;

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
  [ '~/temp/medulla_align_test_global_', num2str(config.stack.case_ids(1)), ...
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
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Superpixel - multi stage
%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(length(config.superpixel)==1)
%   config.superpixel = repmat(config.superpixel, [1 1]);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.1 Superpixel - agglomerative mean boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel(1).method = 'seeded_watershed';
% () thresholds on the boundary field
config.superpixel(1).seeded_watershed.seed_filter = 't0.3_neg_mito0\_e5\_a40_d2_e2_a20';
% Whether to save the segmentation overlaps as TIF files for viewing.
config.superpixel(1).save_segmentation_overlay = true;
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel(1).use_mitochondria = false;
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
config.superpixel(1).watershed_filter_version = 'db_nwo700_cbf3_LDA7_mf7_ps1\_0.25_abf_neg';
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel(1).filter_version = 'db_nwo700_cbf3_LDA7_mf7_ps1\_0.25_abf_neg';
% () whether to print messages
config.superpixel(1).is_verbose = is_verbose;
% () whether to display intermediate results
config.superpixel(1).is_verbose_figures = is_verbose_figures;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7 Choose a superpixel segmentation parameter for further processing
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two options:
% --- (a) One parameter for all images ---
% % Superpixel method
config.superpixel_choose.choice.method = '#(SP_CHOICE)';
% % use superpixel segmentation parameters from previous step
% % E.g., for grayscale-ladder it would be of the form
% % .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria m<threshold_d<amount of erosion>>  
config.superpixel_choose.choice.seg_suffix = '#(SP_CHOICE)';
% --- (b) For each image choose a parameter from a set ---
% % Use a GUI to choose between a set of parameters for each image.
% % --- (c) Set a parameter for a batch of sections ---
% % Useful if blocks of sections have different parameters.
config.superpixel_choose.set_choice.method = 'seeded_watershed';
config.superpixel_choose.set_choice.seg_suffix = ...
  [ '.sws', ...
    '.db_nwo700_cbf3_LDA7_mf7_ps1_0.25_abf_neg_f', ...
    '.db_nwo700_cbf3_LDA7_mf7_ps1_0.25_abf_neg'];

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.1 Superpixel - prune-classify
%%%%%%%%%%%%%%%%%%%%%%%%%%
if(length(config.superpixel_2_seg)==1)
  config.superpixel_2_seg = repmat(config.superpixel_2_seg, [2 1]);
end
% (a) Algorithm to be used - must match with the function called
config.superpixel_2_seg(1).method = 'segment_and';
% () Superpixel method
config.superpixel_2_seg(1).superpixel_method = '';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
config.superpixel_2_seg(1).superpixel_suffix = ''; 
% Which segment maps to merge into the superpixel map
% ** Segmentation method
config.superpixel_2_seg(1).segment_and.to_merge_segment_maps(1).method = 'prune_classify';
% ** Segmentation suffix
config.superpixel_2_seg(1).segment_and.to_merge_segment_maps(1).seg_suffix = ...
  [ '.apc_sp_region_max_3_0.35', ...
    '.gs_l_sp_T0.25_L300', ...
    '.gs_amdb_sp_T0.35_0.34_b5', ...
    '.gs_amb_T0.25_0.24_b0', ...
    '.db_nwo700_cbf3_LDA7_mf7_ps1_0.25_abf_neg_f', ...
    '.db_nwo700_cbf3_LDA7_mf7_ps1_0.25_abf_neg', ...
    '.db_nwo700_cbf3_LDA7_mf7_ps1_0.25_abf_neg_f', ...
    '.db_nwo700_cbf3_LDA7_mf7_ps1_0.25_abf_neg_f', ...
    '.db_nwo700_cbf3_LDA7_mf7_ps1_0.25_abf_neg_f'];
% ** Suffix to be appended to the file name
config.superpixel_2_seg(1).segment_and.to_merge_segment_maps(1).save_suffix = ...
  '.pc_0.35';
% Threshold on the area of overlap between two segments for them to be
% merged (pixel units).
config.superpixel_2_seg(1).segment_and.overlap_area_threshold = 150;
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
config.superpixel_2_seg(1).filter_version = 'db_nwo700_cbf3_LDA7_mf7_ps1\_0.25_abf_neg';
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
config.superpixel_2_seg(2).superpixel_method = 'segment_and';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
config.superpixel_2_seg(2).superpixel_suffix = ...
  [ '.sa_sp_150.pc_0.35', ...
    config.superpixel_choose.choice.seg_suffix, ...
    '.db_nwo700_cbf3_LDA7_mf7_ps1_0.25_abf_neg_f'];
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
config.superpixel_2_seg(2).filter_version = 'db_nwo700_cbf3_LDA7_mf7_ps1\_0.25_abf_neg';
% () whether to print messages
config.superpixel_2_seg(2).is_verbose = is_verbose;
% () whether to display intermediate results
config.superpixel_2_seg(2).is_verbose_figures = is_verbose_figures;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9 Choose a segmentation parameter for further processing
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two options:
% --- (a) One parameter for all images ---
config.segmentation_choose.choice.method = '#(SEG_CHOICE)';
config.segmentation_choose.choice.seg_suffix = '#(SEG_CHOICE)';
% --- (b) For each image choose a parameter from a set ---
% % Use a GUI to choose between a set of parameters for each image.
% --- (c) Set a parameter for a batch of sections ---
% Useful if blocks of sections have different parameters.
config.segmentation_choose.set_choice.method = 'segment_and';
config.segmentation_choose.set_choice.seg_suffix = ...
  [ '.sa_sp_150.pc_0.35', ...
    config.superpixel_choose.choice.seg_suffix, ...
    '.db_nwo700_cbf3_LDA7_mf7_ps1_0.25_abf_neg_f'];

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
config.linkage.feature.filter_version = 'db_nwo700_cbf2_LDA7_mf7_ps1\_0.25_abf_neg';
% Feature suffix for the file name
config.linkage.feature.suffix = ['.', num2str(config.linkage.feature.bin_size), ...
  '.', regexprep(config.linkage.feature.filter_version, ...
  '(?<char>[^\\])\', '$<char>')];
% Use a modified version of segmentation for the linkage - specify prefix
% and suffix to the segmentation parameter string for retrieving the files
config.linkage.modified_segment.method = 'shrink_levelset';
config.linkage.modified_segment.prefix = '.skl_sp_mul5_i25';
config.linkage.modified_segment.suffix = '.db_nwo700_cbf3_LDA7_mf7_ps1_0.25_abf_neg_f';
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
% config.proofreader.roi.xmin = 6511;
% config.proofreader.roi.ymin = 7311;
% config.proofreader.roi.xmax = 12510;
% config.proofreader.roi.ymax = 12310;
% config.proofreader.roi.xmin = 6911;
% config.proofreader.roi.ymin = 6011;
% config.proofreader.roi.xmax = 10410;
% config.proofreader.roi.ymax = 9510;
config.proofreader.roi.xmin = 8411;
config.proofreader.roi.ymin = 7511;
config.proofreader.roi.xmax = 10410;
config.proofreader.roi.ymax = 9510;
config.proofreader.export.is_verbose_figures = is_verbose_figures;
% If grayscale and superpixel maps are to be rendered externally. In this
% case, skip modules 810.1, 810.15, 810.2, 810.3, 810.35. When 810.4 is
% called, only the superpixel-to-segment maps are re-indexed.
% config.proofreader.export.render_externally = true;
% Annotations
% config.proofreader.annotations.method = 'proofread_result';
% config.proofreader.annotations.suffix = '.prfd_annot.concatenate_161_1460';

return
end
