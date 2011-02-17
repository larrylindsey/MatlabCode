function mitochondria_sandbox()

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
reconstruction_param.mitochondria.train.window_sizes = [15 25 35];
reconstruction_param.mitochondria.train.intensity_bins = 0:0.1:1;

reconstruction_param.mitochondria.apply.dir = 'mitochondria/';
reconstruction_param.mitochondria.apply.save_suffix = '.mitochondria_det_conf';

% mitochondria_manual_annotation(reconstruction_param);

% mitochondria_collect_training_samples_intensity_histograms(reconstruction_param);

% mitochondria_train_boosted_intensity_hist_detector(reconstruction_param);

mitochondria_apply_boosted_intensity_hist_detector(reconstruction_param);

return

end