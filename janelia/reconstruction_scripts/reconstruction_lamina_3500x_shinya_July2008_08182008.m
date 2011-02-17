function reconstruction_lamina_3500x_shinya_July2008_08182008(module_id)
% shiv_HPF_FS_S3_R1_001_008_segmentation_08142008(module_id)
% 2D segmentation of stack medulla_HPF_2500x_S3-R1.001_008
% based on optimizations performed on ministack.
% Algo: grayscale ladder for superpixels, followed by Agglo. mean boundary
% + Agglo. median boundary + Agglo. boundary vs interior Fisher Discriminant.
%
% module_id
%   0     SIFT alignment
%   1     mitochondria_manual_annotation(config);
%   2     mitochondria_collect_training_samples_intensity_histograms(config);
%   3     mitochondria_train_boosted_intensity_hist_detector(config);
%   4     mitochondria_apply_boosted_intensity_hist_detector(config);
%   5     segment_2D(config);
%   6     superpixel_2_segment(config);
%   7     linkage_3D_train(config);
%   8     linkage_3D_gen_linkage_graph(config);
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
% v1  05072008  code taken from takemuras/shinya_OsO4_8_10_04292008.m
% v2  06042008  code taken from takemuras/shinya_HPF_FS_000_002_05192008.m
% v3  07152008  code taken from takemuras/shinya_HPF_FS_000_008_06182008.m
% v4  07292008  code taken from takemuras/shinya_HPF_FS_005_007_07152008.m
% v5  07292008  code taken from takemuras/shinya_HPF_FS_S3_R1_001_008_07292008.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. Constants - Not to be changed by users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config = get_basic_config();

config.DEBUG = true;

% Force images to read from Marta's w.s.'s hard disk for faster processing
config.stack.dir = '/opt/Marta/';
config.reconstruction.root_dir = '/opt/Marta/reconstructions/';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% II. Configuration parameters - to be set by user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Stack parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) stack name
config.stack.name = 'data'; % 'lamina.OTO.shinya.leginon.3500x.July2008';   
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
config.stack.image_structure = 'Jcartridge_e.xml';
% * slice ids: specify the z-plane numbers.
config.stack.case_ids = 0; % :3;
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
% (b) Consistency threshold for RANSAC to consider aligning 
% a pair of images
config.SIFT.RANSAC.consistency=0.3;
% () number of iterations of RANSAC. Transformation with minimum error
% among all iterations is chosen.
config.SIFT.n_ransac_iter = 10;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) name for the reconstruction
config.reconstruction.name = 'J_09052008';

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
config.mitochondria.train.case_ids = 1:4;
% For constructing the feature vector for mitochondria detection.
% () type of features collected for each pixel
config.mitochondria.feature.type = 'heq_intensity_hist';
% * version name for the feature
config.mitochondria.feature.version = '_v6';
% * The patch size is (2xwindows_size+1) X (2xwindows_size+1).
config.mitochondria.feature.window_sizes = [4 7 10 13 18 23];
% * In case of intensity histograms, their bins are specified in intensity_bins
config.mitochondria.feature.intensity_bins = 0:0.1:1;
% * The spatial extent of the Gabor filter (gabor_size+1) X (gabor_size+1).
config.mitochondria.feature.gabor.size = 20;
% * The number theta orientations - response is max for all theta.
config.mitochondria.feature.gabor.n_theta = 8;
% * The spatial frequencies for the Gabor filters
config.mitochondria.feature.gabor.freqs = [0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4];
% () Classifier
% * type
config.mitochondria.model.type = 'boost';
% * version
config.mitochondria.model.version = '_b';
% * Number of iterations of boosting
config.mitochondria.model.n_iteration = 30;
% * Tree depth
config.mitochondria.model.tree_depth = 3;
% () Suffix for saving mitochondria detection confidences
config.mitochondria.apply.save_suffix = ['.mitochondria_det_conf.hih', ...
  config.mitochondria.model.version, config.mitochondria.feature.version];
% () Mitochondria detection evaluation
% * version to evaluate
config.mitochondria.evaluate.suffix = config.mitochondria.apply.save_suffix;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Vesicle
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.vesicle.apply.dir = '';
config.vesicle.apply.save_suffix = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Superpixel
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel.method = 'grayscale_ladder';
% (b) thresholds on the boundary field
config.superpixel.f_thresholds = 0.25;
% (c) minimum area of a segment - below this they merged with neighbors
config.superpixel.area_thresholds = 100;
% (d) Whether mitochondria are to be used to suppress false boundaries
config.superpixel.use_mitochondria =false;
% (e) The confidence threshold to be applied on the mitochondria
config.superpixel.mitochondria_confidence_threshold = 0.0;
% (f) The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel.mitochondria_erosion = 5;
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
config.superpixel.filter_version = 'v_mf2';
% () whether to display intermediate results
config.superpixel.is_verbose = true;


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
config.superpixel_2_seg(1).superpixel_method = 'grayscale_ladder';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria m<threshold_d<amount of erosion>>
config.superpixel_2_seg(1).superpixel_suffix = '.gs_l_T0.25_L100.v_mf2'; 
% (c) thresholds on the boundary field
config.superpixel_2_seg(1).f_threshold_seq = 0.26:0.01:0.6;
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
config.superpixel_2_seg(1).mitochondria_erosion = 10;
% () Whether to use vesicle detection to suppress false boundaries
config.superpixel_2_seg(1).use_vesicle = false;
% () The threshold for the vesicle detection confidence
config.superpixel_2_seg(1).vesicle_threshold = 0;
% () Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.1
config.superpixel_2_seg(1).vesicle_dilation = 0;
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel_2_seg(1).filter_version = 'v_heq_mf2';
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
config.superpixel_2_seg(2).superpixel_suffix = '.gs_amb_sp_T0.55_0.54_b0.gs_l_T0.25_L100.v_mf2.v_heq_mf2';
% '.gs_l_T0.47_L100.v_heq_mf3_e2_m0_d2_a500'; 
% (c) thresholds on the boundary field
config.superpixel_2_seg(2).f_threshold_seq = 0.50:0.01:0.62; 
% 0.48:0.01:0.9;
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
config.superpixel_2_seg(2).mitochondria_erosion = 10;
% () Whether to use vesicle detection to suppress false boundaries
config.superpixel_2_seg(2).use_vesicle = false;
% () The threshold for the vesicle detection confidence
config.superpixel_2_seg(2).vesicle_threshold = 0;
% () Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.
config.superpixel_2_seg(2).vesicle_dilation = 0;
% () Image filtering before performing segmentation (affects boundary
% detection)
config.superpixel_2_seg(2).filter_version = 'v_heq_mf2';
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
  '.gs_amdb_sp_T0.6_0.59_b0.gs_amb_sp_T0.55_0.54_b0.gs_l_T0.25_L100.v_mf2.v_heq_mf2.v_heq_mf2'; % '.BEL_l_T0.4_L200.5_b_v_mf5_d2'
% (c) thresholds on the boundary field
config.superpixel_2_seg(3).f_threshold_seq = 0.1:0.01:0.8; % :0.02:0.7; % 0.002:0.07;
% () minimum boundary length (pixel units) for two segments to be merged
% should increase with increasing f_thresholds. E.g., 0 for f_thresholds
% <=0.65 and 20 otherwise
config.superpixel_2_seg(3).length_threshold_seq = ...
  10*(config.superpixel_2_seg(3).f_threshold_seq>0);
% () Whether mitochondria are to be excluded from image statistics
config.superpixel_2_seg(3).exclude_mitochondria_stat = false;
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel_2_seg(3).use_mitochondria = false;
% () The confidence threshold to be applied on the mitochondria
config.superpixel_2_seg(3).mitochondria_confidence_threshold = 0.0;
% () The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel_2_seg(3).mitochondria_erosion = 10;
% () Whether to use vesicle detection to suppress false boundaries
config.superpixel_2_seg(3).use_vesicle = false;
% () The threshold for the vesicle detection confidence
config.superpixel_2_seg(3).vesicle_threshold = 0;
% () Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.
config.superpixel_2_seg(3).vesicle_dilation = 0;
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel_2_seg(3).filter_version = 'v_mf2';
% () whether to display intermediate results
config.superpixel_2_seg(3).is_verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.B Stitch together the segmentations from different tiles
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.align_tile.is_verbose = true;
config.align_tile.superpixel_method = 'grayscale_ladder';
config.align_tile.segment_suffix = ...
    ['.gs_abif_sp_T0.35_0.34_b10', config.superpixel_2_seg(end).superpixel_suffix, '.v_mf2'];
config.align_tile.method = 'boundary_bipartite_match';
config.align_tile.delta_margin = 7;
config.align_tile.min_n_votes = 15;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. 3D linkage graph - training
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) 3D annotation file to be used for training
config.linkage.train.manual_annotation_file = ...
 '/groups/chklovskii/chklovskiilab/electron_microscopy_data/lamina.OTO.3500x.zhiyuan.May232008.3_5/manual_annotations/marta.06032008.wcat.mat';  
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
['.gs_abif_sp_T0.35_0.34_b10', config.superpixel_2_seg(end).superpixel_suffix, '.v_mf2'];
% (j) Verbose
config.linkage.apply.is_verbose = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8. Final reconstruction volume
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
