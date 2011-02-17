function copy_proofread_roi_to_output_dir(config)
% copy_proofread_roi_to_output_dir(config)
% copy proofread ROI to data_to_be_proofread directory for records.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

if(config.is_verbose)
  fprintf('START: copy_proofread_roi_to_output_dir\n');
end
proofreader_roi = config.proofreader.roi; %#ok<NASGU>
proofreader_roi_file_name = ...
  [get_global_stitching_param_dir(config), config.proofreader.roi.file_name];
save2(proofreader_roi_file_name, 'proofreader_roi');
proofreader_roi_file_name = ...
  [get_to_be_proofread_dir(config), 'global_stitching_parameters/', ...
  config.proofreader.roi.file_name];
save2(proofreader_roi_file_name, 'proofreader_roi');

if(config.is_verbose)
  fprintf('STOP: copy_proofread_roi_to_output_dir\n');
end
return
end
