function output_to_raveler_stitched_grayscale(al, config)
% output_to_raveler_stitched_grayscale(al, config)
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

stack_config = config.stack;
raveler_config = config.proofreader.Raveler;

save_dir = [get_to_be_proofread_dir(config), raveler_config.grayscale_dir];
check_for_dir(save_dir);

for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  file_name = [save_dir, 'image.', raveler_config.grayscale_version_name, '.', ...
    num2str(case_id, '%05d'), '.png'];
  imwrite(uint8(al{i}), file_name, 'BitDepth', 8);
end

return
end
