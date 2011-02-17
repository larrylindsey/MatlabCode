function reconstruction_sandbox()
% reconstruction_sandbox()
% Testbed for the reconstruction pipeline
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stack parameters
reconstruction_param.name = '04122008.boost_02262008';
reconstruction_param.dir = ['./data_to_be_proofread/Marta_aligned/', reconstruction_param.name, '/'];
reconstruction_param.image_dir = '/groups/chklovskii/chklovskiilab/Shiv/Marta_aligned/';
reconstruction_param.image_prefix = 'a%04d';
reconstruction_param.image_suffix = '.tif';
reconstruction_param.case_ids = 1:4; % 20
reconstruction_param.roi.xmin = 2134;
reconstruction_param.roi.xmax = 4412;
reconstruction_param.roi.ymin = 1817;
reconstruction_param.roi.ymax = 4095;


% mitochondria
reconstruction_param.mitochondria.train.dir = 'mitochondria/training/';
reconstruction_param.mitochondria.train.image_prefix = 'a%04d_crop';
reconstruction_param.mitochondria.train.image_suffix = '.tif';
reconstruction_param.mitochondria.train.save_suffix = '_mitochondria.mat';
reconstruction_param.mitochondria.train.case_ids = [1, 2];
reconstruction_param.mitochondria.train.model_name = 'mitochondria_boost_detect.mat';
reconstruction_param.mitochondria.train.window_sizes = [15 25 35];
reconstruction_param.mitochondria.train.intensity_bins = 0:0.1:1;

reconstruction_param.mitochondria.apply.dir = 'mitochondria/';
reconstruction_param.mitochondria.apply.save_suffix = '.mitochondria_det_conf';


% vesicle
reconstruction_param.vesicle.apply.dir = '';
reconstruction_param.vesicle.apply.save_suffix = '';


% superpixel
reconstruction_param.superpixel_param.dir = '2D_segmentation_results/grayscale_ladder/';
reconstruction_param.superpixel_param.use_mitochondria = true;
reconstruction_param.superpixel_param.mitochondria_confidence_threshold = 0.0;
reconstruction_param.superpixel_param.use_vesicle = false;
reconstruction_param.superpixel_param.f_thresholds = 0.49;
reconstruction_param.superpixel_param.area_thresholds = 400;


% superpixel to segment
reconstruction_param.superpixel_2_seg_param.dir = '2D_segmentation_results/grayscale_ladder/';
reconstruction_param.superpixel_2_seg_param.superpixel_suffix = '.gs_l_T0.49_L400.v_heq_mf1_m0'; % use superpixel segmentation parameters
reconstruction_param.superpixel_2_seg_param.use_mitochondria = false;
reconstruction_param.superpixel_2_seg_param.mitochondria_confidence_threshold = 0.0;
reconstruction_param.superpixel_2_seg_param.use_vesicle = false;
reconstruction_param.superpixel_2_seg_param.f_thresholds = 0.55;
reconstruction_param.superpixel_2_seg_param.area_thresholds = 400;


% 3D linkage
reconstruction_param.linkage.train.dir = 'linkage_3D/';
reconstruction_param.linkage.train.manual_annotation_file = ...
  '/groups/chklovskii/chklovskiilab/electron_microscopy_data/fly.C155-Elav_UAS-CD2-HRP.5000x.zhiyuan.Feb052008/manual_annotations/annotation_2D_seg_3D_link.Shiv_02252008.mat';
reconstruction_param.linkage.train.intensity_bins = 0:0.1:1.0;
reconstruction_param.linkage.train.save_suffix = '.04112008';

reconstruction_param.linkage.apply.dir = 'linkage_3D/';
reconstruction_param.linkage.apply.model_suffix = '.04112008';
reconstruction_param.linkage.apply.segmentation_suffix = ... % use segmentation parameters
  ['.gs_l_sp_T0.55_L400', reconstruction_param.superpixel_2_seg_param.superpixel_suffix, '.v_heq_mf1'];


% save_reconstruction_param_to_xml(reconstruction_param);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segmentation and Linkage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mitochondria_manual_annotation(reconstruction_param);

% mitochondria_collect_training_samples_intensity_histograms(reconstruction_param);

% mitochondria_train_boosted_intensity_hist_detector(reconstruction_param);

% mitochondria_apply_boosted_intensity_hist_detector(reconstruction_param);

% reconstruction_param.segmentation_2D = reconstruction_param.superpixel_param;
% segment_2D_grayscale_ladder(reconstruction_param);
% reconstruction_param = rmfield(reconstruction_param, 'segmentation_2D');

reconstruction_param.segmentation_2D = reconstruction_param.superpixel_2_seg_param;
superpixel_2_segment_grayscale_ladder(reconstruction_param);
reconstruction_param = rmfield(reconstruction_param, 'segmentation_2D');

% linkage_3D_train_intensity_pair_boost(reconstruction_param);

% linkage_3D_gen_linkage_gph_intensity_pair_boost(reconstruction_param);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the final data structures to be given to the proofreding gui.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate_al(reconstruction_param);

% generate_cat(reconstruction_param);

% generate_superpixel_2_seg_map(reconstruction_param);

% process_superpixel_2_2Dseg_linkage_graph(reconstruction_param);

% save_reconstruction_info(reconstruction_param);

return

end