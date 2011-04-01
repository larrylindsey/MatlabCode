function config = mitochondria_2007KHLateral_s09_1_91(module_id)
% config = mitochondria_2007KHLateral_s09_1_91(module_id)
% Mitochondria detection on 2007KH Lateral S09 images slices 1 through 91.
%
% module_id
%   1     mitochondria_manual_annotation(config);
%   2     mitochondria_collect_training_samples_intensity_histograms(config);
%   3.2   Train and apply mitochondria detector on the training data using
%         a Leave-one-out protocol. Calls:
%         mitochondria_train_boosted_intensity_hist_detector(config);
%         mitochondria_apply_boosted_intensity_hist_detector(config);
%   4     mitochondria_apply_boosted_intensity_hist_detector(config);
%   4.2   mitochondria_evaluate_detector(config);
%   4.3   mitochondria_evaluate_leave_one_out(config);
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
% v1  05012008  code borrowed from reconstruction_ex_Alex_volume.m
% v2  05132008  code borrowed from mitochondria_OsO4_1_medulla_ihg.m
%

if(nnz(ismember(module_id, [1 2 3 3.2 4 4.2 4.3]))<length(module_id))
  error('Invalid module id/s');
end;

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
config.stack.name = '2007KHLateral.s09.cropped';
% (b) image name:
% If the files are image_0001.tif, image_0002.tif, etc., then
% the image_prefix would be "image_%04d" and the image_suffix would be
% ".tif".
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
config.stack.image_prefix = 's09%03d.crop';  
config.stack.image_suffix = '.tif';
% (c) slice ids: If the images are a001.tif, ..., a010.tif then the
% case_ids are 1,2,...,10.
config.stack.case_ids = 2:90;
% (d) ROI: Region Of Interest for the stack. If the entire slice image is
% to be processed then assign empty ([]). Otherwise get_roi(config);
config.stack.roi = [];
% config.stack.roi = get_roi(config);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Reconstruction 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) name for the reconstruction
config.reconstruction.name = 'mitochondria';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Mitochondria
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Training image name:
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
config.mitochondria.train.image_prefix = 'm.%03d';
config.mitochondria.train.image_suffix = '.tif';
% (b) image ids: : If the images are a006.tif, ..., a010.tif then the
% case_ids are 6,...,10.
config.mitochondria.train.case_ids = 15:10:45;
% For constructing the feature vector for mitochondria detection.
% () type of features collected for each pixel
config.mitochondria.feature.type = 'heq_intensity_hist_3s2';
% * version name for the feature
config.mitochondria.feature.version = '_v5';
% * The patch size is (2xwindows_size+1) X (2xwindows_size+1).
config.mitochondria.feature.window_sizes = [4 7 10 15 20 25 30];
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
config.mitochondria.model.version = '_a3';
% * Number of iterations of boosting
config.mitochondria.model.n_iteration = 100;
% * Tree depth
config.mitochondria.model.tree_depth = 1;
% () Suffix for saving mitochondria detection confidences
config.mitochondria.apply.save_suffix = ['.mitochondria_det_conf.hih3s2', ...
  config.mitochondria.model.version, config.mitochondria.feature.version];
% () Mitochondria detection evaluation
% * version to evaluate
config.mitochondria.evaluate.suffix = config.mitochondria.apply.save_suffix;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. Pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pipeline(module_id, config);

return

end
