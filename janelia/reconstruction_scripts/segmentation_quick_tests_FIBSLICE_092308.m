function segmentation_quick_tests_FIBSLICE_092308(module_id)
% segmentation_quick_tests_FIBSLICE_092308(module_id)
% 2D segmentation of quick_tests images
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
% v6  08182008  code taken from shiv_HPF_FS_S3_R1_001_008_segmentation_08142008.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. Constants - Not to be changed by users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config = get_basic_config();

config.stack.dir = '~/research/em_reconstruction_pipeline/data_em_images/';

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
config.stack.name = 'quick_tests';
% (b) image name:
% If the files are image_0001.tif, image_0002.tif, etc., then
% the image_prefix would be "image_%04d" and the image_suffix would be
% ".tif".
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
config.stack.image_prefix = 'FIBSLICE%04d';  
config.stack.image_suffix = '.tif';
% (c) slice ids: If the images are a001.tif, ..., a010.tif then the
% case_ids are 1,2,...,10.
config.stack.case_ids = 5:8;
% (d) ROI: Region Of Interest for the stack. If the entire slice image is
% to be processed then assign empty ([]). Otherwise get_roi(config);
config.stack.roi = [];
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
config.stack.align.xf_file_name = '';
config.stack.align.xg_file_name = '';
config.stack.align.margin = 200;
config.stack.align.roi = [];
% (f) Stack resolution parameters for area and boundary histograms. These
% are to be used if the models having been trained on images of resolution
% different from the current one. This might improve the results,
% especially during boot-strapping.
config.stack.length_factor = 3;
config.stack.area_factor = config.stack.length_factor.^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) name for the reconstruction
config.reconstruction.name = '09232008';

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
config.superpixel.method = 'grayscale_AglMeanB';
% () thresholds on the boundary field
config.superpixel.f_threshold_seq = 0.1:0.01:0.9;
% () minimum boundary length (pixel units) for two segments to be merged
% should increase with increasing f_thresholds. E.g., 0 for f_thresholds
% <=0.65 and 20 otherwise
config.superpixel.length_threshold_seq = ...
  25*(config.superpixel.f_threshold_seq>0.65);
% () Whether mitochondria are to be used to suppress false boundaries
config.superpixel.use_mitochondria = false;
% () The confidence threshold to be applied on the mitochondria
config.superpixel.mitochondria_confidence_threshold = 0; 
% () The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel.mitochondria_erosion = 5; % 20
% * Minimum of a thresholded connected component for it to be considered
% mitochondria
config.superpixel.mitochondria_min_area = 40;
% () Whether to use vesicle detection to suppress false boundaries
config.superpixel.use_vesicle = false;
% () The threshold for the vesicle detection confidence
config.superpixel.vesicle_threshold = 0;
% () Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.
config.superpixel.vesicle_dilation = 0;
% Image filtering for boundary detection to be used for watershed
% boundary = 1 - filtered_image
config.superpixel.watershed_filter_version = 'v_heq_o1';
% () Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel.filter_version = 'v_heq';
% () whether to display intermediate results
config.superpixel.is_verbose = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6 Multistage superpixel to segment computation
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.superpixel_2_seg = repmat(config.superpixel_2_seg, [1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.1 Superpixel to segment - agglomerative mean boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel_2_seg(1).method = 'grayscale_ladder';
% () Superpixel method
config.superpixel_2_seg(1).superpixel_method = 'grayscale_AglMeanB';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
config.superpixel_2_seg(1).superpixel_suffix = '.gs_amb_T0.81_0.8_b25.v_heq_o1.v_heq'; 
% () thresholds on the boundary field
config.superpixel_2_seg(1).f_thresholds = 0.1;
% () minimum area of a segment - below this they merged with neighbors
config.superpixel_2_seg(1).area_thresholds =50;
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
config.superpixel_2_seg(1).filter_version = 'v_heq_mf3';
% () whether to display intermediate results
config.superpixel_2_seg(1).is_verbose = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. 3D linkage graph - training
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) 3D annotation file to be used for training
config.linkage.train.manual_annotation_file = ...
  '/groups/chklovskii/chklovskiilab/electron_microscopy_data/medulla_HPF_2500x_S3-R1.005_007/manual_annotations/proofread.shinya.07282008.wcat.mat';
% (b) Version name for the linkage model
config.linkage.train.save_suffix = '.08122008';
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
config.linkage.apply.model_suffix = '.08122008';
% (i) The segmenation to be used for linkage
% use the segmentation parameters from previous step
% E.g., for grayscale-ladder
% .gs_l_T<f_threshold>_L<area_threshold>.[superpixel_suffix].<filter_version>_<mitochondria>
config.linkage.apply.segmentation_suffix = ... 
  ['.gs_l_sp_T0.1_L50', config.superpixel_2_seg(1).superpixel_suffix, '.v_heq_mf3']; % '.v_mf5_d1'
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
