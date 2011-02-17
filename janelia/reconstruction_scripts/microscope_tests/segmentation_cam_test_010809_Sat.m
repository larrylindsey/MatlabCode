function segmentation_cam_test_010809_Sat(module_id)
% segmentation_cam_test_010809_Sat(module_id)
% 2D segmentation of stack carbon_thickness_tests_120808_2.2nm for
% evaluating EM camera. Sat vs. 0.25Sat
% Algo: grayscale ladder for superpixels, followed by Agglo. mean boundary
% + Agglo. median boundary + Agglo. boundary vs interior Fisher Discriminant.
%
% module_id
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
% v6  10012008  to test scintillator
% v7  12082008  to test carbon thickness
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. Constants - Not to be changed by users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config = get_basic_config();
global config_global
% Whether to use hashing when generating file names during storage and
% retrieval. This reduces the filename lengths.
config_global.storage.file_name_is_hashed__ = false;

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
config.stack.name = 'cam_test_010809_Spirit_Exposure';   
% (b) image name:
% If the files are image_0001.tif, image_0002.tif, etc., then
% the image_prefix would be "image_%04d" and the image_suffix would be
% ".tif".
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
config.stack.image_prefix = '27932M1.Sat.bin%d.4800x';  
config.stack.image_suffix = '.tif';
% (c) slice ids: If the images are a001.tif, ..., a010.tif then the
% case_ids are 1,2,...,10.
config.stack.case_ids = 2;
% (d) ROI: Region Of Interest for the stack. If the entire slice image is
% to be processed then assign empty ([]). Otherwise get_roi(config);
config.stack.roi = [];
%config.stack.roi.xmin = 200;
%config.stack.roi.xmax = 1800;
%config.stack.roi.ymin = 200;
%config.stack.roi.ymax = 1800;
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
config.stack.align.is_prealigned = true;
% * If not prealigned then name the .xf file specifying the alignment.
% Assumed to be located in the stack directory.
config.stack.align.xf_file_name = '??.xf';
config.stack.align.xg_file_name = '??.xg';
config.stack.align.margin = 200;
config.stack.align.roi = [];
% (f) Stack resolution parameters for area and boundary histograms. These
% are to be used if the models having been trained on images of resolution
% different from the current one. This might improve the results,
% especially during boot-strapping.
config.stack.length_factor = 1;
config.stack.area_factor = config.stack.length_factor.^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) name for the reconstruction
config.reconstruction.name = '08122008';

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
config.mitochondria.train.case_ids = 5:7;
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
% 5. Superpixel
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel.method = 'grayscale_ladder';
% (b) thresholds on the boundary field
config.superpixel.f_thresholds = 0.3:0.02:0.5;
% (c) minimum area of a segment - below this they merged with neighbors
config.superpixel.area_thresholds = 75;
% (d) Whether mitochondria are to be used to suppress false boundaries
config.superpixel.use_mitochondria = false;
% (e) The confidence threshold to be applied on the mitochondria
config.superpixel.mitochondria_confidence_threshold = 0.0;
% (f) The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel.mitochondria_erosion = 2;
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
config.superpixel.filter_version = 'v_heq_mf2';
% () whether to display intermediate results
config.superpixel.is_verbose = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6 Multistage superpixel to segment computation
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.superpixel_2_seg = repmat(config.superpixel_2_seg, [2 1]);
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
% config.superpixel_2_seg(1).superpixel_suffix = ...
%   '.gs_l_T0.36_L75.v_heq_mf2'; % 0.25Sat.4800xbin2
config.superpixel_2_seg(1).superpixel_suffix = ...
  '.gs_l_T0.44_L75.v_heq_mf2'; % Sat.4800xbin2
% (c) thresholds on the boundary field
config.superpixel_2_seg(1).f_threshold_seq = 0.44:0.01:0.9;
% (c) thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel_2_seg(1).save_f_thresholds = 0.6:0.02:0.9;
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
config.superpixel_2_seg(1).filter_version = 'v_heq_mf2';
% () whether to display intermediate results
config.superpixel_2_seg(1).is_verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.2 Superpixel to segment - agglomerative boundary vs. interior values
% Fisher Linear Discriminant
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel_2_seg(2).method = 'grayscale_AglBIFL';
% () Superpixel method
config.superpixel_2_seg(2).superpixel_method = 'grayscale_AglMeanB';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
% For 1900x
% config.superpixel_2_seg(2).superpixel_suffix = ...
%   '.gs_amb_sp_T0.78_0.77_b25.gs_l_T0.36_L75.v_heq_mf2.v_heq_mf2'; % 0.25Sat.4800xbin2
config.superpixel_2_seg(2).superpixel_suffix = ...
  '.gs_amb_sp_T0.78_0.77_b25.gs_l_T0.44_L75.v_heq_mf2.v_heq_mf2'; % Sat.4800xbin2
% (c) thresholds on the boundary field
config.superpixel_2_seg(2).f_threshold_seq = 0.01:0.01:0.3; % :0.02:0.7; % 0.002:0.07;
% (c) thresholds on the boundary field for which segmentations should be
% saved. Useful for reducing disk space.
config.superpixel_2_seg(2).save_f_thresholds = 0.3;
% () minimum boundary length (pixel units) for two segments to be merged
% should increase with increasing f_thresholds. E.g., 0 for f_thresholds
% <=0.65 and 20 otherwise
config.superpixel_2_seg(2).length_threshold_seq = ...
  10*(config.superpixel_2_seg(2).f_threshold_seq>0);
% Maximum area threshold: segments with area above this threshold are not
% considered for merging. Intuition: Segments with very large areas are
% likely to have complex structures with multimodal intensity
% distributions. Merging them is more likely to result in errors.
config.superpixel_2_seg(2).max_area_threshold = 5000;
% () Whether mitochondria are to be excluded from image statistics
config.superpixel_2_seg(2).exclude_mitochondria_stat = false;
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
config.superpixel_2_seg(2).filter_version = 'v_mf3';
% () whether to display intermediate results
config.superpixel_2_seg(2).is_verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. 3D linkage graph - training
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) 3D annotation file to be used for training
config.linkage.train.manual_annotation_file = ...
  '/groups/chklovskii/chklovskiilab/electron_microscopy_data/fly_larva.no_ua.5000xbinx2/manual_annotations/wayne.06032008.wcat.mat';
% (b) Version name for the linkage model
config.linkage.train.save_suffix = '.09092008';
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
config.linkage.apply.model_suffix = '.09092008';
% (i) The segmenation to be used for linkage
% use the segmentation parameters from previous step
% E.g., for grayscale-ladder
% .gs_l_T<f_threshold>_L<area_threshold>.[superpixel_suffix].<filter_version>_<mitochondria>
config.linkage.apply.segmentation_suffix = ... 
  ['.gs_abif_sp_T0.04_0.03_b10', config.superpixel_2_seg(end).superpixel_suffix, '.v_mf3'];
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
