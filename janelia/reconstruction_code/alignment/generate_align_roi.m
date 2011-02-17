function generate_align_roi(config)
% generate_align_roi(config)

align_roi = get_roi(config); %#ok<NASGU>
align_roi_file_name = [get_global_stitching_param_dir(config), ...
  config.stack.align.roi_file_name];
save2(align_roi_file_name, 'align_roi');
align_roi_file_name = [get_to_be_proofread_dir(config), 'global_stitching_parameters/'...
  config.stack.align.roi_file_name];
save2(align_roi_file_name, 'align_roi');
return
end
