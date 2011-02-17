function output_to_raveler_stitched_superpixel_maps(cat, config)
% output_to_raveler_stitched_superpixel_maps(cat, config)
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

stack_config = config.stack;
raveler_config = config.proofreader.Raveler;

save_dir = [get_to_be_proofread_dir(config), raveler_config.superpixel_dir];
check_for_dir(save_dir);

for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  fprintf('plane: %d\n', case_id);
  file_name = [save_dir, 'superpixel_map.', ...
    raveler_config.superpixel_version_name, '.', num2str(case_id, '%05d'), ...
    '.png'];
  imwrite(uint16(cat{i}), file_name, 'BitDepth', 16);
end

return
end
