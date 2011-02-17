function segmentation_sandbox()

reconstruction_param.dir = '/groups/chklovskii/chklovskiilab/Shiv/Marta_aligned/';
reconstruction_param.image_prefix = 'a%04d';
reconstruction_param.image_suffix = '.tif';
reconstruction_param.case_ids = [3, 4];
reconstruction_param.roi.xmin = 2134;
reconstruction_param.roi.xmax = 4412;
reconstruction_param.roi.ymin = 1817;
reconstruction_param.roi.ymax = 4095;


reconstruction_param.mitochondria.train.dir = 'mitochondria/training/';
reconstruction_param.mitochondria.train.image_prefix = 'a%04d_crop';
reconstruction_param.mitochondria.train.image_suffix = '.tif';
reconstruction_param.mitochondria.train.save_suffix = '_mitochondria.mat';
reconstruction_param.mitochondria.train.case_ids = [1, 2];
reconstruction_param.mitochondria.train.model_name = 'mitochondria_boost_detect.mat';

reconstruction_param.mitochondria.apply.dir = 'mitochondria/';
reconstruction_param.mitochondria.apply.save_suffix = '.mitochondria_det_conf';


reconstruction_param.vesicle.apply.dir = '';
reconstruction_param.vesicle.apply.save_suffix = '';


reconstruction_param.superpixel_param.use_mitochondria = true;
reconstruction_param.superpixel_param.use_vesicle = false;
reconstruction_param.superpixel_param.mitochondria_confidence_threshold = 0.0;
reconstruction_param.superpixel_param.f_thresholds = 0.49;
reconstruction_param.superpixel_param.area_thresholds = 400;
reconstruction_param.superpixel_param.dir = '2D_segmentation_results/grayscale_ladder/';

reconstruction_param.segmentation_2D = reconstruction_param.superpixel_param;


reconstruction_param.superpixel_2_seg_param.use_mitochondria = false;
reconstruction_param.superpixel_2_seg_param.use_vesicle = false;
reconstruction_param.superpixel_2_seg_param.mitochondria_confidence_threshold = 0.0;
reconstruction_param.superpixel_2_seg_param.f_thresholds = 0.55;
reconstruction_param.superpixel_2_seg_param.area_thresholds = 400;
reconstruction_param.superpixel_2_seg_param.dir = '2D_segmentation_results/grayscale_ladder/';
reconstruction_param.superpixel_2_seg_param.superpixel_suffix = '_gs_l_T0.49_L400.v_heq_mf1_m0'; % use superpixel segmentation parameters

reconstruction_param.segmentation_2D = reconstruction_param.superpixel_2_seg_param;


reconstruction_param.linkage.train.intensity_bins = 0:0.1:1.0;
reconstruction_param.linkage.train.manual_annotation_file = ...
  '/groups/chklovskii/chklovskiilab/electron_microscopy_data/fly.C155-Elav_UAS-CD2-HRP.5000x.zhiyuan.Feb052008/manual_annotations/annotation_2D_seg_3D_link.Shiv_02252008.mat';
reconstruction_param.linkage.train.dir = 'linkage_3D/';
reconstruction_param.linkage.train.save_suffix = '.04112008';

reconstruction_param.linkage.apply.segmentation_suffix = ... % use segmentation parameters
  ['_gs_l_sp_T0.55_L400', reconstruction_param.superpixel_2_seg_param.superpixel_suffix, '.v_heq_mf1'];


linkage_3D_train_intensity_pair_boost(reconstruction_param);


return

end