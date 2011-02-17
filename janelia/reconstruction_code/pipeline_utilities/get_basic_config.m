function config = get_basic_config()
% get_basic_config()
% Gives basic config parameters that should usually not be changed by
% users.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. Constants - Not to be changed by users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global config_global

config.is_verbose = true;
config.is_verbose_figures = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Global configuration parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Whether to use hashing when generating file names during storage and
% retrieval. This reduces the filename lengths.
config_global.storage.file_name_is_hashed__ = true;
% Important: Don't change the following values
config_global.storage.hash.string_index_cycle = 4;
config_global.storage.hash.multiplication_factor = (sqrt(5)-1)/2;
config_global.storage.hash.hash_range = 2^20;
% Check to make hash function is compatibile
if(~check_hash_compatibility())
  error('Hash function is not compatible. Saved file names will not be determined.');
end
% During storage, if the file names are to be hashed then part of the file
% name string is hashed. The following define the prefix and suffix of the
% file name that are kept intact for legibility.
% Important: Don't change the following values
config_global.storage.hash.leave_as_is_prefix_length = 50;
config_global.storage.hash.min_length_to_hash = 20;
config_global.storage.hash.leave_as_is_suffix_length = 4;
% For deformable mesh based alignment, the valid mesh region labels
% (transform ids) start from TRANSFORMATION_ID_OFFSET. Values 0 to
% TRANSFORMATION_ID_OFFSET are reserved for special regions such as folds,
% etc. 
config_global.TRANSFORMATION_ID_OFFSET = 10;
% For generating name prefixs from a 2-tuple of strings. Useful for
% generating file names for objects defined over tuples, e.g., transforms
% between pairs of images. 
config_global.name_prefix_from_string_tuple.leave_as_is_prefix_length = 6;
% Whether running in debug mode. Behavior:
%   print debug messages
config_global.DEBUG = true;
% Whether to execute routines within MATLAB or generate stand alone batch
% scripts that % can be run directly in Linux. Set job.is_stand_along to
% true for generating Linux stand-alone scripts.
if(~isfield(config_global, 'job') || ...
    ~isfield(config_global.job, 'is_stand_alone'))
  config_global.job.is_stand_alone = false; % default
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Stack parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Directory within which the stack directories are stored
config.stack.dir = '/groups/chklovskii/chklovskiilab/electron_microscopy_data/';
config.stack.image_suffix = '.tif';
config.stack.align.xf_save_name = 'xf.mat';
config.stack.align.xg_save_name = 'xg.mat';
config.stack.align.roi_file_name = 'align_roi.mat';
% (b) Stack resolution parameters for area and boundary histograms. These
% are to be used if the models having been trained on images of resolution
% different from the current one. This might improve the results,
% especially during boot-strapping. Default values:
config.stack.length_factor = 1;
config.stack.area_factor = 1;

config.region.registry_dir = 'region_registry/';

config.job.dir = 'jobs/';
config.job.log.deformable_mesh_in_plane = 'deformable_mesh_in_plane.log';
config.job.log.deformable_mesh_inter_plane = 'deformable_mesh_inter_plane.log';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Directory within which the recontruction output files are stored
prev_dir = pwd2();
curr_dir = prev_dir;
found_reconstruction_dir = false;
while((strcmp(curr_dir, '/')~=1) && ~found_reconstruction_dir)
  cd('..');
  curr_dir = pwd2();
  if(exist('reconstructions', 'dir')==7)
    found_reconstruction_dir = true;
  end
end
if(~found_reconstruction_dir)
  error('Not in correct directory. Please go within reconstructions/scripts/ and run again.');
end
config.reconstruction.root_dir = [pwd2, '/'];
addpath([config.reconstruction.root_dir, 'scripts/']);
cd(prev_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Mitochondria
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Directory for storing mitochondria specific files
config.mitochondria.dir = 'mitochondria/';
% (b) Directory for training information
config.mitochondria.train.dir = 'training/';
% (c) Suffix for save mitochondria manual annotations
config.mitochondria.train.annot_suffix = '.mitochondria_annot';
% (d) Mitochondria model suffix
config.mitochondria.train.model_prefix = 'mitochondria_model';
% (e) Suffix for saving the mitochondria detection confidences
config.mitochondria.apply.save_suffix = '.mitochondria_det_conf';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Vesicle
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.vesicle.dir = 'vesicles/';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Superpixel
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.superpixel.dir = '2D_segmentation_results/';
config.superpixel.save_segmentation_overlay = true;

config.BEL.dir = 'BEL/';

config.LDA.dir = 'LDA/';
config.LDA.train.dir = 'models/';

config.precompute.filtered_image.dir = 'filtered_images/';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Superpixel to segment
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.segmentation_2D.dir = '2D_segmentation_results/';
config.segmentation_2D.save_segmentation_overlay = true;
config.superpixel_2_seg.dir = '2D_segmentation_results/';
config.superpixel_2_seg.save_segmentation_overlay = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. 3D linkage graph
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.linkage.dir = 'linkage_3D/';
config.linkage.evaluate.segmentation_suffix = '.eval_k';
config.linkage.evaluate.temporary_dir = '/tmp/';
config.linkage.evaluate.save_dir = 'eval_results/';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8. Final reconstruction volume
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.final_volume.dir = 'seg_volume/';

config.to_be_proofread.root_dir = ...
  '/groups/chklovskii/chklovskiilab/em_reconstruction/data_to_be_proofread/';
config.proofreading_session.root_dir = ...
  '/groups/chklovskii/chklovskiilab/em_reconstruction/proofreading_sessions/';

config.coarse_al.file_name = 'al_coarse.tif';
config.coarse_al.scale = 4;
config.coarse_al.dir = 'al_coarse/';
config.coarse_al.sectionwise_file_name_format = 'al_coarse.%07d.tif';
config.coarse_al.save_sections_separately = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9. Evaluate reconstruction volume
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.evaluate_volume.annotations_dir = 'manual_annotations/';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Evaluate segment boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.eval_segment_boundary.dir = 'bp_match_results/';
config.eval_segment_boundary.annotations_dir = 'manual_annotations/';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11. Training segmentation algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.train_segmentation.dir = 'train_segmentation/';
config.train_segmentation.annotations_dir = 'manual_annotations/';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12. Superpixel parameter choice
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.superpixel_choose.dir = 'superpixel_param_choices/';
config.superpixel_choose.input_keys = {'a', 'p', 'b', 'd', 'f', 'g', 'h'};
config.superpixel_choose.plot_colors = {'r', 'y', 'b', 'g'};
config.superpixel_choose.display_top_half = true;
config.superpixel_choose.save_suffix = '.superpixel_param';
config.superpixel_choose.is_image_histeq = true;


%%%
% Desired decision tree - some default trees defined for convenience
%%%
% This specifies the manner in which choices are to be presented to the
% user.
% Type of decision, e.g., 2 for binary is two choices are to be presented.
config.decision_tree_predef_binary_4_choice.type = 2;
% The decision tree consists of a tree of nodes. Each node is a struct with
% fields:
%   .param_id_0   seg. params for "Left" choice
%   .param_id_1   seg. params for "Right choice
%   .noide_id_0   tree node id to go to if user chooses "Left"
%   .noide_id_1   tree node id to go to if user chooses "Right"
%
% The first node is always the root node - the first choice presented to
% the user.
config.decision_tree_predef_binary_4_choice.node(1).param_id(1) = 2;
config.decision_tree_predef_binary_4_choice.node(1).param_id(2) = 3;
config.decision_tree_predef_binary_4_choice.node(1).node_id(1) = 2;
config.decision_tree_predef_binary_4_choice.node(1).node_id(2) = 3;

config.decision_tree_predef_binary_4_choice.node(2).param_id(1) = 1;
config.decision_tree_predef_binary_4_choice.node(2).param_id(2) = 2;
config.decision_tree_predef_binary_4_choice.node(2).node_id(1) = -1;
config.decision_tree_predef_binary_4_choice.node(2).node_id(2) = -2;

config.decision_tree_predef_binary_4_choice.node(3).param_id(1) = 3;
config.decision_tree_predef_binary_4_choice.node(3).param_id(2) = 4;
config.decision_tree_predef_binary_4_choice.node(3).node_id(1) = -3;
config.decision_tree_predef_binary_4_choice.node(3).node_id(2) = -4;

% Type of decision, e.g., 2 for binary is two choices are to be presented.
config.decision_tree_predef_binary_4_choice_bias_4.type = 2;
% The decision tree consists of a tree of nodes. Each node is a struct with
% fields:
%   .param_id_0   seg. params for "Left" choice
%   .param_id_1   seg. params for "Right choice
%   .noide_id_0   tree node id to go to if user chooses "Left"
%   .noide_id_1   tree node id to go to if user chooses "Right"
%
% The first node is always the root node - the first choice presented to
% the user.
config.decision_tree_predef_binary_4_choice_bias_4.node(1).param_id(1) = 1;
config.decision_tree_predef_binary_4_choice_bias_4.node(1).param_id(2) = 4;
config.decision_tree_predef_binary_4_choice_bias_4.node(1).node_id(1) = 2;
config.decision_tree_predef_binary_4_choice_bias_4.node(1).node_id(2) = -4;

config.decision_tree_predef_binary_4_choice_bias_4.node(2).param_id(1) = 2;
config.decision_tree_predef_binary_4_choice_bias_4.node(2).param_id(2) = 3;
config.decision_tree_predef_binary_4_choice_bias_4.node(2).node_id(1) = 3;
config.decision_tree_predef_binary_4_choice_bias_4.node(2).node_id(2) = -3;

config.decision_tree_predef_binary_4_choice_bias_4.node(3).param_id(1) = 1;
config.decision_tree_predef_binary_4_choice_bias_4.node(3).param_id(2) = 2;
config.decision_tree_predef_binary_4_choice_bias_4.node(3).node_id(1) = -1;
config.decision_tree_predef_binary_4_choice_bias_4.node(3).node_id(2) = -2;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 13. Segemntation parameter choice
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.segmentation_choose.dir = 'segmentation_param_choices/';
config.segmentation_choose.input_keys = {'a', 'p', 'b', 'd', 'f', 'g', 'h'};
config.segmentation_choose.plot_colors = {'r', 'y', 'b', 'g'};
config.segmentation_choose.display_top_half = true;
config.segmentation_choose.save_suffix = '.segmentation_param';
config.segmentation_choose.is_image_histeq = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 14. Normalized cuts - images segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.normalized_cuts.dir = 'normalized_cuts/';
config.normalized_cuts.filter_name_prefix = 'ncut_';
% maximum side of image to be given to ncuts. Larger images are downsampled
config.normalized_cuts.max_image_side = 400;
% number of segments, approximately equal to number of eigen vectors to be
% computed by normalized cuts 
config.normalized_cuts.n_segments = 60;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 15. Alignment of stack - multiple tiles in a section, pairs of tiles in
% adjacent section, and multiple tiles multiple sections.
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.align_segmentation.dir = 'align_segment_map/';

config.deformable_mesh.dir = 'deformable_mesh/';

config.norm_x_corr.dir = 'norm_x_corr/';
config.align.precompute.norm_cross_corr_in_plane.scale = 5;
config.align.precompute.norm_cross_corr_inter_plane.scale = 5;

config.SIFT.dir = 'alignmentSIFT/';
% For SIFT RANSAC
% matches with deviation under current model within this factor of current
% inliers' mean deviation to be added to inliers
config.align.linkage_align.SIFT.RANSAC.inlier_max_relative_deviation = 1;
% fraction of best matches under current model to be considered inliers
config.align.linkage_align.SIFT.RANSAC.best_inlier_carryover_frac = 0.15;
% allowed discrepancy between consistent matches
config.align.linkage_align.SIFT.RANSAC.max_discrepancy_for_consistency = 16;
% minimal number of RANSAC iterations
config.align.linkage_align.SIFT.RANSAC.min_n_iter = 10;
% amount of noise added to vote count of consistency. This induces a
% certain amount of randomness to the ranking, and thus the sampling 
config.align.linkage_align.SIFT.RANSAC.rank_perturb_noise = 2;
% Format to save/load global patchwise affine transforms [default '.mat']
config.align.global_align.storage_format = 'mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 16. Folds
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.fold.dir = 'fold_masks/';
% Suffix for storing fold masks
config.fold.save_suffix = '.fold_mask';
% Minimum area of a component in the thresholded image for it to be
% considered as a fold.
config.precompute.fold.min_fold_area = 20000; % default
% Whether to save the fold masks as TIF files that can be used by
% non-MATLAB programs
config.fold.save_as_tif = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 17. Output datastructures to trakEM
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directory for storing the patch masks
config.trakEM.output_patch_masks_transforms.dir = config.fold.dir;
% File specifying the preamble to the java-script. Copied verbatim
config.trakEM.output_patch_masks_transforms.preamble_file = ...
  'trakEM_output_patch_masks_transforms.preamble.txt';
% File specifying the postscript to the java-script. Copied verbatim
config.trakEM.output_patch_masks_transforms.postscript_file = ...
  'trakEM_output_patch_masks_transforms.postscript.txt';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 18. Output datastructures to proofreader
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.proofreader.method = 'matlab_gui'; % default

% Whether to simply generate an identity ROI for the stack. If set to
% false, convex hulls are computed for all sections and a vertical overlap
% mask is computed.
config.proofreader.align_roi.is_identity = true;

config.proofreader.export = [];
% If the superpixel maps are being rendered externally then the remapping
% is saved in these text files.
config.proofreader.export.superpixel_remapping_file = 'map-sp.%d.txt';

config.proofreader.Raveler.grayscale_dir = 'grayscale_maps/';
config.proofreader.Raveler.grayscale_version_name = 'v_a';
config.proofreader.Raveler.superpixel_dir = 'superpixel_maps/';
config.proofreader.Raveler.superpixel_version_name = 'v_a';
config.proofreader.Raveler.superpixel_to_segment_file_name = ...
  'superpixel_to_segment_map.txt';
config.proofreader.Raveler.segment_to_body_map_file_name = ...
  'segment_to_body_map.txt';

config.proofreader.import.dir = 'proofread_result/';
config.proofreader.import.Raveler.superpixel_dir = 'sp_maps/';
config.proofreader.import.Raveler.superpixel_prefix = 'sp_map.%d';
config.proofreader.import.Raveler.superpixel_suffix = '.png';
config.proofreader.import.Raveler.superpixel_to_segment_file_name = ...
  'superpixel-to-segment-map.txt';
config.proofreader.import.Raveler.segment_to_body_map_file_name = ...
  'segment-to-body-map.txt';
config.proofreader.roi.file_name = 'proofreader_roi.mat';
config.proofreader.import.Raveler.annotation_file_name = 'annotations.pickle';

config.annotations.dir = 'annotations/';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 19. 3D segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.image_stack.dir = 'image_stacks/';

config.segmentation_3D.dir = '3D_segmentation_results/';

config.initial_segmentation_3D = config.segmentation_3D;
config.superpixel_3D = config.segmentation_3D;

config.combine_segment_3D.dir = 'combined_seg_3D/';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 20. Notification email
%%%%%%%%%%%%%%%%%%%%%%%%%%
config_global.notify_email.gmailUser = 'em.reconstruct';
config_global.notify_email.gmailPassword = 'janelia.em.reconstruct';
config_global.notify_email.recipient = 'vitaladevunis@janelia.hhmi.org';

return
end
