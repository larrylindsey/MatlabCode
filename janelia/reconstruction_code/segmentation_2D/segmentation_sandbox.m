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

% superpixel_param.use_mitochondria = true;
% superpixel_param.use_vesicle = false;
% superpixel_param.mitochondria_confidence_threshold = 0.0;
% superpixel_param.f_thresholds = 0.49;
% superpixel_param.area_thresholds = 400;
% superpixel_param.dir = '2D_segmentation_results/grayscale_ladder/';
% 
% reconstruction_param.segmentation_2D = superpixel_param;
% 
% segment_2D_grayscale_ladder(reconstruction_param);

superpixel_2_seg_param.use_mitochondria = false;
superpixel_2_seg_param.use_vesicle = false;
superpixel_2_seg_param.mitochondria_confidence_threshold = 0.0;
superpixel_2_seg_param.f_thresholds = 0.55;
superpixel_2_seg_param.area_thresholds = 400;
superpixel_2_seg_param.dir = '2D_segmentation_results/grayscale_ladder/';
superpixel_2_seg_param.superpixel_suffix = '_gs_l_T0.49_L400.v_heq_mf1_m0';

reconstruction_param.segmentation_2D = superpixel_2_seg_param;

superpixel_2_segment_grayscale_ladder(reconstruction_param);

return

end