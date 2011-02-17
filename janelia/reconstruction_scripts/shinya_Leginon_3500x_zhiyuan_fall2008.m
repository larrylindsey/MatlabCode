function shinya_Leginon_3500x_zhiyuan_fall2008(module_id)
% test_downsampled_seg_medulla_HPF_leginon(module_id)
% Test downsampled reconstruction on medulla stack
% chklovskii/medulla/medulla.HPF.Leginon.3500x.zhiyuan.fall2008
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
% v3  10142008  modified for medulla
% v4  10302008  adding ABIF within superpixel to reduce oversegmentation in
%                 dark neurites
% v5  12032008  save files in the chklovskii/medulla/ directory
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. Constants - Not to be changed by users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config = get_basic_config();

config.stack.dir = '/groups/chklovskii/medulla/';
% Uncomment to save intermediate reconstruction files in the
% chklovskii/medulla/ directory
config.reconstruction.root_dir = '/groups/chklovskii/medulla/reconstructions/';

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
% config.stack.image_prefix = 'a.%03d';  
% config.stack.image_suffix = '.tif';
% If the images are a001.tif, ..., a010.tif then the
% case_ids would be 1:10.
% config.stack.case_ids = 1:20;
% Option 2: If the images correspond to a trakEM project then
% * specify the trakEM file structure xml file name
config.stack.image_structure = 'crop2.ms.8bit.3500x9x9_0161_0170.xml';
% * slice ids: specify the z-plane numbers.
config.stack.case_ids = 161:170;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.1. SIFT alignment
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Calculate adjacency matrix from trakEM2 tiling
config.SIFT.trakEM2.getAdjacency=true;
% (b) Consistency threshold for RANSAC to consider aligning 
% a pair of images
config.SIFT.RANSAC.consistency=0.3;
% () number of iterations of RANSAC. Transformation with minimum error
% among all iterations is chosen.
config.SIFT.n_ransac_iter = 100;
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
config.SIFT.is_inverted_display = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) name for the reconstruction
config.reconstruction.name = 'test_SIFT_alignment';

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
config.superpixel = repmat(config.superpixel, [3 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.1 Superpixel - agglomerative mean boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel(1).method = 'grayscale_AglMeanB';
% () thresholds on the boundary field
config.superpixel(1).f_threshold_seq = 0.1:0.01:0.6; % 0.1:0.01:0.9;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel(1).save_f_thresholds = [0.45,0.5,0.55,0.58,0.6]; %0.5:0.02:0.6;
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
config.superpixel(1).watershed_filter_version = 'v_heq_o1';
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel(1).filter_version = 'v_heq2';
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
config.superpixel(2).superpixel_suffix = '.gs_amb_T0.58_0.57_b0.v_heq_o1_m0_d5_a40.v_heq2'; %'.gs_amb_T0.6_0.59_b0.v_heq_o1_m0_d5_a40.v_heq'; 
% () thresholds on the boundary field
config.superpixel(2).f_thresholds = 0.2;
% () minimum area of a segment - below this they merged with neighbors
config.superpixel(2).area_thresholds =60;
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
config.superpixel(2).filter_version = 'v_heq2';
% () whether to display intermediate results
config.superpixel(2).is_verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.3 Superpixel - agglomerative boundary vs. interior values
% Fisher Linear Discriminant
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel(3).method = 'grayscale_AglBIFL';
% () Superpixel method
config.superpixel(3).superpixel_method = 'grayscale_ladder';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
config.superpixel(3).superpixel_suffix = ...
  '.gs_l_sp_T0.2_L60.gs_amb_T0.58_0.57_b0.v_heq_o1_m0_d5_a40.v_heq2.v_heq2';
% (c) thresholds on the boundary field
config.superpixel(3).f_threshold_seq = 0.1:0.01:0.14; % :0.02:0.7; % 0.002:0.07;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel(3).save_f_thresholds = 0.14; %0.14:0.02:0.2;
% () minimum boundary length (pixel units) for two segments to be merged
% should increase with increasing f_thresholds. E.g., 0 for f_thresholds
% <=0.65 and 20 otherwise
config.superpixel(3).length_threshold_seq = ...
  2*(config.superpixel(3).f_threshold_seq>0);
% () Maximum area of a segment, above this segments are not considered for
% merging.
config.superpixel(3).max_area_threshold = 40000;
% () Whether mitochondria are to be excluded from image statistics
config.superpixel(3).exclude_mitochondria_stat = true;
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel(3).use_mitochondria = false;
% () The confidence threshold to be applied on the mitochondria
config.superpixel(3).mitochondria_confidence_threshold = 0.0;
% () The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel(3).mitochondria_erosion = 20;
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
config.superpixel(3).filter_version = 'v0';
% () whether to display intermediate results
config.superpixel(3).is_verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6 Multistage superpixel to segment computation
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.superpixel_2_seg = repmat(config.superpixel_2_seg, [3 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.1 Superpixel to segment - agglomerative mean boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel_2_seg(1).method = 'grayscale_AglMeanB';
% () Superpixel method
config.superpixel_2_seg(1).superpixel_method = 'grayscale_AglBIFL';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria m<threshold_d<amount of erosion>>
config.superpixel_2_seg(1).superpixel_suffix = ''; %'.gs_l_sp_T0.2_L60.gs_amb_T0.6_0.59_b0.v_heq_o1_m0_d5_a40.v_heq2.v_heq2_mf2'; 
% (c) thresholds on the boundary field
config.superpixel_2_seg(1).f_threshold_seq = 0.4:0.01:0.64; %0.4:0.01:0.7;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel_2_seg(1).save_f_thresholds = 0.52:0.02:0.64;
% () minimum boundary length (pixel units) for two segments to be merged
% should increase with increasing f_thresholds. E.g., 0 for f_thresholds
% <=0.65 and 20 otherwise
config.superpixel_2_seg(1).length_threshold_seq = ...
  25*(config.superpixel_2_seg(1).f_threshold_seq>0.65);
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
config.superpixel_2_seg(1).filter_version = 'v_heq2';
% () whether to display intermediate results
config.superpixel_2_seg(1).is_verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.2 Superpixel to segment - agglomerative median boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel_2_seg(2).method = 'grayscale_AglMedianB';
% () Superpixel method
config.superpixel_2_seg(2).superpixel_method = 'grayscale_AglMeanB';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria m<threshold_d<amount of erosion>>
config.superpixel_2_seg(2).superpixel_suffix = '.gs_amb_sp_T0.52_0.51_b0.gs_abif_sp_T0.14_0.13_b2.gs_l_sp_T0.2_L60.gs_amb_T0.52_0.51_b0.v_heq_o1_m0_d5_a40.v_heq2.v_heq2.v0.v_heq2';
% '.gs_l_T0.47_L100.v_heq_mf3_e2_m0_d2_a500'; 
% (c) thresholds on the boundary field
config.superpixel_2_seg(2).f_threshold_seq = 0.5:0.01:0.64;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel_2_seg(2).save_f_thresholds = 0.52:0.02:0.64;
% () minimum boundary length (pixel units) for two segments to be merged
% should increase with increasing f_thresholds. E.g., 0 for f_thresholds
% <=0.65 and 20 otherwise
config.superpixel_2_seg(2).length_threshold_seq = ...
  25*(config.superpixel_2_seg(2).f_threshold_seq>0.65);
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
config.superpixel_2_seg(2).filter_version = 'v_heq2';
% () whether to display intermediate results
config.superpixel_2_seg(2).is_verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.3 Superpixel to segment - agglomerative boundary vs. interior values
% Fisher Linear Discriminant
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel_2_seg(3).method = 'grayscale_AglBIFL';
% () Superpixel method
config.superpixel_2_seg(3).superpixel_method = 'grayscale_AglMedianB';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
config.superpixel_2_seg(3).superpixel_suffix = ...
  '.gs_amdb_sp_T0.52_0.51_b0.gs_amb_sp_T0.52_0.51_b0.gs_abif_sp_T0.14_0.13_b2.gs_l_sp_T0.2_L60.gs_amb_T0.52_0.51_b0.v_heq_o1_m0_d5_a40.v_heq2.v_heq2.v0.v_heq2.v_heq2';
% (c) thresholds on the boundary field
config.superpixel_2_seg(3).f_threshold_seq = 0.1:0.01:0.7; % :0.02:0.7; % 0.002:0.07;
% () thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel_2_seg(3).save_f_thresholds = 0.2:0.02:0.5;
% () minimum boundary length (pixel units) for two segments to be merged
% should increase with increasing f_thresholds. E.g., 0 for f_thresholds
% <=0.65 and 20 otherwise
config.superpixel_2_seg(3).length_threshold_seq = ...
  10*(config.superpixel_2_seg(3).f_threshold_seq>0);
% () Maximum area of a segment, above this segments are not considered for
% merging.
config.superpixel_2_seg(3).max_area_threshold = 40000;
% () Whether mitochondria are to be excluded from image statistics
config.superpixel_2_seg(3).exclude_mitochondria_stat = true;
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel_2_seg(3).use_mitochondria = false;
% () The confidence threshold to be applied on the mitochondria
config.superpixel_2_seg(3).mitochondria_confidence_threshold = 0.0;
% () The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel_2_seg(3).mitochondria_erosion = 20;
% () Whether to use vesicle detection to suppress false boundaries
config.superpixel_2_seg(3).use_vesicle = false;
% () The threshold for the vesicle detection confidence
config.superpixel_2_seg(3).vesicle_threshold = 0;
% () Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.
config.superpixel_2_seg(3).vesicle_dilation = 0;
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel_2_seg(3).filter_version = 'v0';
% () whether to display intermediate results
config.superpixel_2_seg(3).is_verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. 3D linkage graph - training
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
config.linkage.train.save_suffix = '.11032008';
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
config.linkage.apply.model_suffix = '.11032008';
% () The superpixel map to be used for linkage - must match with the
% segmentation being used.
config.linkage.apply.superpixel_method = 'grayscale_ladder';
% use the segmentation parameters from previous step
% E.g., for grayscale-ladder
% .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
config.linkage.apply.superpixel_suffix = config.superpixel_2_seg(1).superpixel_suffix;
% (i) The segmenation to be used for linkage
config.linkage.apply.segmentation_method = 'align_segment_map';
% use the segmentation parameters from previous step
% E.g., for grayscale-ladder
% .gs_l_T<f_threshold>_L<area_threshold>.[superpixel_suffix].<filter_version>_<mitochondria>
config.linkage.apply.segmentation_suffix = ... 
['.gs_abif_sp_T0.22_0.21_b10', config.superpixel_2_seg(end).superpixel_suffix, '.v0'];
% (j) Verbose
config.linkage.apply.is_verbose = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8. Stitch together the segmentations of different tiles within a section
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.align_tile.is_verbose = true;
config.align_tile.segment_suffix = ...
  ['.gs_abif_sp_T0.22_0.21_b10', config.superpixel_2_seg(end).superpixel_suffix, '.v0'];
config.align_tile.method = 'boundary_bipartite_match';
config.align_tile.delta_margin = 7;
config.align_tile.min_n_votes = 15;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9. Final reconstruction volume
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.final_volume.link_thresholds = -.1:0.1:2;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9. Evaluate reconstruction volume
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.evaluate_volume.groundtruth_file = 'manual_annotation.1_5.Alex.mat';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. Reconstruction Routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pipeline(module_id, config);

return

end
