function test_SIFT_alignment_JC_102208(module_id)
% test_SIFT_alignment_JC_102208(module_id)
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

config.stack.dir = '/groups/chklovskii/chklovskiilab/JC_data/alignment_test/';

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
config.stack.name = 'RGBs';   
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
config.stack.image_structure = 'RGBChannels1-300DSPS_JFImaging.xml';
% * slice ids: specify the z-plane numbers.
config.stack.case_ids = 0:20;
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
% 1.1. SIFT alignment
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Calculate adjacency matrix from trakEM2 tiling
config.SIFT.trakEM2.getAdjacency=true;
% () Filter the image before applying the SIFT detector
config.SIFT.filter_version = 'mul100';
% (b) Consistency threshold for RANSAC to consider aligning 
% a pair of images
config.SIFT.RANSAC.consistency=0.3;
% () number of iterations of RANSAC. Transformation with minimum error
% among all iterations is chosen.
config.SIFT.n_ransac_iter = 50;
% () number of initial sample correspondences used for computing an initial
% estimate of model. Should be greater 2 for affine transformations.
% config.SIFT.n_initial_samples = 6;
% (c) Target downsampled image size for SIFT
config.SIFT.downsampling=512;
% (d) Calculate SIFT landmarks and their descriptors
config.SIFT.makeFrames=true;
% (e) Calculate matches between SIFT landmarks
config.SIFT.makeMatches=true;
% (f) Calculate global affine transforms for each tile
config.SIFT.makeTransforms=true;
% (g) Save SIFT landmarks, matches and transforms
config.SIFT.save=true;
% () Whether to print intermediate messages
config.SIFT.is_verbose = true;
% () Minimum number of matches between two tiles for their matches to be
% included in the section alignment. For pairs of tiles with marginal
% overlap, the matches are dominated by false correspondences. Therefore,
% the purity of union of matches is improved by eliminating small sets.
config.SIFT.min_n_matches = 50;
% () Whether to show the images inverted for easier assessment of alignment
config.SIFT.is_inverted_display = false;
% () Whether to save the displayed images and if so the file name with path
% as desired. Has to be TIF format.
config.SIFT.save_as_stack_file_name = '~/temp/JC_0_20_110408.tif';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. Reconstruction Routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pipeline(module_id, config);

return

end
