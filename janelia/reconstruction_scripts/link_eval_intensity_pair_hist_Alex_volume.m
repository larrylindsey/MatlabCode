function link_eval_intensity_pair_hist_Alex_volume(module_id)
% link_eval_intensity_pair_hist_Alex_volume(module_id)
% Evaluate linkage algorithm - intensity pair histogram - on Alex volume.
% Pairs of segments are used for training and evaluation in a n-fold
% protocol.
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
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. Constants - Not to be changed by users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config = get_basic_config();

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
config.stack.name = 'recon_ex';   
% (b) image name:
% If the files are image_0001.tif, image_0002.tif, etc., then
% the image_prefix would be "image_%04d" and the image_suffix would be
% ".tif".
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
config.stack.image_prefix = 'a%03d';  
config.stack.image_suffix = '.tif';
% (c) slice ids: If the images are a001.tif, ..., a010.tif then the
% case_ids are 1,2,...,10.
config.stack.case_ids = 1:5;
% (d) ROI: Region Of Interest for the stack. If the entire slice image is
% to be processed then assign empty ([]). Otherwise get_roi(config);
config.stack.roi = [];
% config.stack.roi = get_roi(config);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) name for the reconstruction
config.reconstruction.name = '04142008';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Mitochondria
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Training image name:
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
config.mitochondria.train.image_prefix = 'a%03d';
config.mitochondria.train.image_suffix = '.tif';
% (b) image ids: : If the images are a006.tif, ..., a010.tif then the
% case_ids are 6,...,10.
config.mitochondria.train.case_ids = 6:10;
% (c) For constructing the feature vector for mitochondria detection. In
% case intensity histograms, the patch size is (2xwindows_size+1) X (2xwindows_size+1).
% The intensity histogram's bins are specified in intensity_bins
config.mitochondria.feature.type = 'heq_intensity_hist';
config.mitochondria.feature.window_sizes = [15 25 35];
config.mitochondria.feature.intensity_bins = 0:0.1:1;
% (d) Classifier type
config.mitochondria.model.type = 'boost';
% (e) Number of iterations of boosting
config.mitochondria.model.n_iteration = 30;
% (f) Tree depth
config.mitochondria.model.tree_depth = 1;

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
config.superpixel.f_thresholds = 0.56;
% (c) minimum area of a segment - below this they merged with neighbors
config.superpixel.area_thresholds = 600;
% (d) Whether mitochondria are to be used to suppress false boundaries
config.superpixel.use_mitochondria = true;
% (e) The confidence threshold to be applied on the mitochondria
config.superpixel.mitochondria_confidence_threshold = 0.0;
% (f) The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel.mitochondria_erosion = -1;
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
config.superpixel.filter_version = 'v_heq_mf10';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Superpixel to segment
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel_2_seg.method = 'grayscale_ladder';
% (b) Superpixel version to be used as basis
% % use superpixel segmentation parameters from previous step
% For grayscale-ladder .gs_l_T<f_threshold>_L<area_threshold>.<filter_version>_<mitochondria>
config.superpixel_2_seg.superpixel_suffix = '.gs_l_T0.56_L600.v_heq_mf10_m0'; 
% (c) thresholds on the boundary field
config.superpixel_2_seg.use_mitochondria = false;
% (d) minimum area of a segment - below this they merged with neighbors
config.superpixel_2_seg.mitochondria_confidence_threshold = 0.0;
% (e) Whether mitochondria are to be used to suppress false boundaries
config.superpixel_2_seg.mitochondria_erosion = 20;
% (f) The confidence threshold to be applied on the mitochondria
config.superpixel_2_seg.use_vesicle = false;
% (g) The thresholded mitochondria mask is eroded to ensure that correct
% boundaries are not obliterated
config.superpixel_2_seg.vesicle_threshold = 0;
% (h) Whether to use vesicle detection to suppress false boundaries
config.superpixel_2_seg.vesicle_dilation = 0;
% (i) The threshold for the vesicle detection confidence
config.superpixel_2_seg.f_thresholds = 0.58;
% (j) Radius of the vesicles. Boundaries within this distance from detected
% vesicle centers are obliterated.
config.superpixel_2_seg.area_thresholds = 600;
% (k) Image filtering before performing segmentation (affects boundary
% detection).
config.superpixel_2_seg.filter_version = 'v_heq_mf10';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. 3D linkage graph - training
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) 3D annotation file to be used for training
config.linkage.train.manual_annotation_file = ...
  '/groups/chklovskii/chklovskiilab/Shiv/Alex_volume/linkage_3D/training/proofread_6_10.Alex.mat';
% (b) Version name for the linkage model
config.linkage.train.save_suffix = '.04242008d';
% (c) Type for feature to be used for linkage
config.linkage.feature.type = 'intensity_pair_hist';
% (d) For intensity pair hist., the intensity bins in the histograms
config.linkage.feature.intensity_bins = 0:0.1:1.0;
% (e) Type of classifier
config.linkage.model.type = 'boost';
% (f) Number of iterations of boosting
config.linkage.model.n_iteration = 15;
% (g) Tree depth
config.linkage.model.tree_depth = 2;
% () Number of iterations of evaluation. In each iteration, the bodies are
% randomly permuted and then evaluated using k-fold protocol.
config.linkage.evaluate.n_iteration = 1;
% () k-fold evaluation of linkage algorithm. (k-1)/k fraction of data is used
% for training, the remaining used for testing.
config.linkage.evaluate.k_fold = 5;
% () Range of linkage thresholds for the precision recall cure
config.linkage.evaluate.link_thresholds = -1:0.2:3;
% () Whether to dump crops of slice images of linkage pairs 
% link_dump.<z, label_1, label_2, segment_1, segment_2>.
config.linkage.evaluate.is_link_dumped = true;

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