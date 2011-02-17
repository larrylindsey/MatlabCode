function save_reconstruction_info(config)
% save_reconstruction_info(config)
% save the reconstruction information - to be used when starting the
% proofreading session
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
%

fprintf('Saving reconstruction info. .. ');

linkage_config = config.linkage;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
train_config = linkage_config.train;
apply_config = linkage_config.apply;

i.config = config;
i.al_file = get_storage_file_name([get_to_be_proofread_dir(config), 'al.mat']);
i.cat_file = get_storage_file_name([get_to_be_proofread_dir(config), ...
  'cat', config.superpixel_choose.choice.seg_suffix, '.mat']);
i.superpixel_2_seg_map_file = get_storage_file_name(...
  [get_to_be_proofread_dir(config), 'superpixel_2_seg_map', ...
  config.segmentation_choose.choice.seg_suffix, '.mat']);
i.links_3D_file = get_storage_file_name([get_to_be_proofread_dir(config), ...
  'links_3D',  '.',  model_config.type, '.', feature_config.type, ...
  apply_config.model_suffix, config.segmentation_choose.choice.seg_suffix, '.mat']);
i.link_model_file = get_storage_file_name([get_to_be_proofread_dir(config), ...
  'link_model', '.', model_config.type, '.', ...
  feature_config.type, train_config.save_suffix, '.mat']);
  
save2([get_to_be_proofread_dir(config), config.reconstruction.name, '.to_be_proofread.mat'], '-STRUCT', 'i');

fprintf('done.\n');

