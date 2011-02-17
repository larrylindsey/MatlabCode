function reconstruction_davi_data(module_id)
% reconstruction_davi_data(module_id)
% reconstruction davi data
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
% v1  06052008
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
config.stack.name = 'davi_data_06052008';   
% (b) image name:
% If the files are image_0001.tif, image_0002.tif, etc., then
% the image_prefix would be "image_%04d" and the image_suffix would be
% ".tif".
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
% config.stack.image_prefix = '03_cardona_grid1_section01_col4_row3_im%d_1_00x_minification.crop1'; 
% config.stack.case_ids = 18;
% config.stack.image_prefix = '080516_sect63_g5A_ms7_5L_450x900_col14_row23_im%d_1_00x_minification'; 
% config.stack.case_ids = 720;
% config.stack.image_prefix = '080516_sect63_g5A_ms7_5L_450x900_col14_row7_im%d_1_00x_minification'; 
% config.stack.case_ids = 736;
config.stack.image_prefix = '080516_sect63_g5A_ms7_5L_450x900_col14_row53_im%d_1_00x_minification'; 
config.stack.case_ids = 690;
config.stack.image_suffix = '.tif';
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
% 4. Vesicle
%%%%%%%%%%%%%%%%%%%%%%%%%%
config.vesicle.apply.dir = '';
config.vesicle.apply.save_suffix = '';

% %%%%%%%%%%%%%%%%%%%%%%%%%%
% % 5. Superpixel
% '03_cardona_grid1_section01_col4_row3_im%d_1_00x_minification.crop1'
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% % (a) Algorithm to be used - must match with the function called
% config.superpixel.method = 'grayscale_ladder';
% % (b) thresholds on the boundary field
% config.superpixel.f_thresholds = 0.44:0.01:0.5;
% % (c) minimum area of a segment - below this they merged with neighbors
% config.superpixel.area_thresholds = 200;
% % (d) Whether mitochondria are to be used to suppress false boundaries
% config.superpixel.use_mitochondria = false;
% % (e) The confidence threshold to be applied on the mitochondria
% config.superpixel.mitochondria_confidence_threshold = 0.0;
% % (f) The thresholded mitochondria mask is eroded to ensure that correct
% % boundaries are not obliterated
% config.superpixel.mitochondria_erosion = -1;
% % * Minimum of a thresholded connected component for it to be considered
% % mitochondria
% config.superpixel.mitochondria_min_area = 500;
% % (g) Whether to use vesicle detection to suppress false boundaries
% config.superpixel.use_vesicle = false;
% % (h) The threshold for the vesicle detection confidence
% config.superpixel.vesicle_threshold = 0;
% % (i) Radius of the vesicles. Boundaries within this distance from detected
% % vesicle centers are obliterated.
% config.superpixel.vesicle_dilation = 0;
% % (j) Image filtering before performing segmentation (affects boundary
% % detection).
% config.superpixel.filter_version = 'v_heq_mf5_e2';

% %%%%%%%%%%%%%%%%%%%%%%%%%%
% % 5. Superpixel
% % '080516_sect63_g5A_ms7_5L_450x900_col14_row23_im720_1_00x_minification' 
% % '080516_sect63_g5A_ms7_5L_450x900_col14_row7_im736_1_00x_minification'
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% % (a) Algorithm to be used - must match with the function called
% config.superpixel.method = 'grayscale_ladder';
% % (b) thresholds on the boundary field
% config.superpixel.f_thresholds = [0.43 0.47 0.5];
% % (c) minimum area of a segment - below this they merged with neighbors
% config.superpixel.area_thresholds = 800;
% % (d) Whether mitochondria are to be used to suppress false boundaries
% config.superpixel.use_mitochondria = false;
% % (e) The confidence threshold to be applied on the mitochondria
% config.superpixel.mitochondria_confidence_threshold = 0.0;
% % (f) The thresholded mitochondria mask is eroded to ensure that correct
% % boundaries are not obliterated
% config.superpixel.mitochondria_erosion = -1;
% % * Minimum of a thresholded connected component for it to be considered
% % mitochondria
% config.superpixel.mitochondria_min_area = 500;
% % (g) Whether to use vesicle detection to suppress false boundaries
% config.superpixel.use_vesicle = false;
% % (h) The threshold for the vesicle detection confidence
% config.superpixel.vesicle_threshold = 0;
% % (i) Radius of the vesicles. Boundaries within this distance from detected
% % vesicle centers are obliterated.
% config.superpixel.vesicle_dilation = 0;
% % (j) Image filtering before performing segmentation (affects boundary
% % detection).
% config.superpixel.filter_version = 'v_heq_mf5_e3';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Superpixel
% '080516_sect63_g5A_ms7_5L_450x900_col14_row53_im690_1_00x_minification'
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Algorithm to be used - must match with the function called
config.superpixel.method = 'grayscale_ladder';
% (b) thresholds on the boundary field
config.superpixel.f_thresholds = 0.43;
% (c) minimum area of a segment - below this they merged with neighbors
config.superpixel.area_thresholds = 3600;
% (d) Whether mitochondria are to be used to suppress false boundaries
config.superpixel.use_mitochondria = false;
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
config.superpixel.filter_version = 'v_heq_mf5_d6';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. Reconstruction Routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pipeline(module_id, config);

return

end