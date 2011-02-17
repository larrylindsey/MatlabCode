function choose_superpixel_config_medulla_HPF_leginon_fall08(module_id)
% choose_superpixel_config_recon_ex(module_id)
% Configuration for choosing superpixel parameters. 
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  11192008  init code
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. Constants - Not to be changed by users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config = get_basic_config();

config.stack.dir = '/groups/chklovskii/medulla/';
% Uncomment to save intermediate reconstruction files in the
% chklovskii/medulla/ directory
% config.reconstruction.root_dir = '/groups/chklovskii/medulla/reconstructions/';

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
config.stack.name = 'medulla.HPF.Leginon.3500x.zhiyuan.fall2008';   
% (b) image name:
% If the files are image_0001.tif, image_0002.tif, etc., then
% the image_prefix would be "image_%04d" and the image_suffix would be
% ".tif".
% If the files are a001.tif, a002.tif, etc., then
% the image_prefix would be "a%03d" and the image_suffix would be
% ".tif".
% Option 2: If the images correspond to a trakEM project then
% * specify the trakEM file structure xml file name
config.stack.image_structure = 'crop2.ms.8bit.3500x9x9_0681_0690.xml';
% * slice ids: specify the z-plane numbers.
config.stack.case_ids = 684:690;
% (d) ROI: Region Of Interest for the stack. If the entire slice image is
% to be processed then assign empty ([]). Otherwise get_roi(config);
config.stack.roi = [];
% config.stack.roi = get_roi(config);
% () Set the scale at which segmentation should be performed. This is kept
% constant for all superpixel and superpixel-to-seg stages. When stitching
% the segmentations, the label maps can be further downsampled if needed,
% but not artificially upsampled. During linkage and proofreading, the
% segmentations are upsampled to original resolution.
config.stack.segmentation_scale = 2;
config.stack.segmentation_scaling_method = 'bilinear';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Superpixel choice parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a list of superpixel parameters. Each element is a struct with
% fields:
%   .method       the segmentation method, e.g., 'grayscale_ladder'.
%   .seg_suffix   parameters for segmentation, similar to those used in
%                   superpixel_suffix in reconstruction scripts.
config.superpixel_choose.param(1).method = 'grayscale_AglBIFL';
config.superpixel_choose.param(1).seg_suffix = ...
  '.gs_abif_sp_T0.18_0.17_b2.gs_l_sp_T0.2_L60.gs_amb_T0.45_0.44_b0.v_heq_o1_m0_d5_a40.v_heq2.v_heq2.v0';

config.superpixel_choose.param(2).method = 'grayscale_AglBIFL';
config.superpixel_choose.param(2).seg_suffix = ...
  '.gs_abif_sp_T0.18_0.17_b2.gs_l_sp_T0.2_L60.gs_amb_T0.45_0.44_b0.v_heq_o1_m0_d5_a40.v_heq2.v_heq2.v0';

config.superpixel_choose.param(3).method = 'grayscale_AglBIFL';
config.superpixel_choose.param(3).seg_suffix = ...
  '.gs_abif_sp_T0.18_0.17_b2.gs_l_sp_T0.2_L60.gs_amb_T0.63_0.62_b0.v_heq_o1_m0_d5_a40.v_heq2.v_heq2.v0';

config.superpixel_choose.param(4).method = 'grayscale_AglBIFL';
config.superpixel_choose.param(4).seg_suffix = ...
  '.gs_abif_sp_T0.18_0.17_b2.gs_l_sp_T0.2_L60.gs_amb_T0.58_0.57_b0.v_heq_o1_m0_d5_a40.v_heq2.v_heq2.v0';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Desired decision tree
%%%%%%%%%%%%%%%%%%%%%%%%%%
% This specifies the manner in which choices are to be presented to the
% user. Two options:
% (a) use one of the predefined decision trees (see get_basic_config.m)
config.superpixel_choose.decision_tree = config.decision_tree_predef_binary_4_choice;
% (b) OR, define a custom decision tree:
% Type of decision, e.g., 2 for binary is two choices are to be presented.
% The decision tree consists of a tree of nodes. Each node is a struct with
% fields:
%   .param_id_0   seg. params for "Left" choice
%   .param_id_1   seg. params for "Right choice
%   .noide_id_0   tree node id to go to if user chooses "Left"
%   .noide_id_1   tree node id to go to if user chooses "Right"
%
% The first node is always the root node - the first choice presented to
% the user.
% See get_basic_config, Sec. 12 for examples.

<<<<<<< .mine
=======
config.superpixel_choose.is_image_histeq = true;
>>>>>>> .r262

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. GUI for choosing superpixel parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pipeline(module_id, config);

return
end
