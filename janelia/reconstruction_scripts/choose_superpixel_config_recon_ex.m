function choose_superpixel_config_recon_ex(module_id)
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
% 1. Superpixel choice parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a list of superpixel parameters. Each element is a struct with
% fields:
%   .method       the segmentation method, e.g., 'grayscale_ladder'.
%   .seg_suffix   parameters for segmentation, similar to those used in
%                   superpixel_suffix in reconstruction scripts.
config.superpixel_choose.param(1).method = 'grayscale_ladder';
config.superpixel_choose.param(1).seg_suffix = '.gs_l_T0.44_L600.v_heq_mf10_m0_d-1_a500';

config.superpixel_choose.param(2).method = 'grayscale_ladder';
config.superpixel_choose.param(2).seg_suffix = '.gs_l_T0.52_L600.v_heq_mf10_m0_d-1_a500';

config.superpixel_choose.param(3).method = 'grayscale_ladder';
config.superpixel_choose.param(3).seg_suffix = '.gs_l_T0.64_L600.v_heq_mf10_m0_d-1_a500';

config.superpixel_choose.param(4).method = 'grayscale_ladder';
config.superpixel_choose.param(4).seg_suffix = '.gs_l_T0.72_L600.v_heq_mf10_m0_d-1_a500';

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

config.superpixel_choose.is_image_histeq = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. GUI for choosing superpixel parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pipeline(module_id, config);

return
end
