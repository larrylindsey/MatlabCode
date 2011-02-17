function output_to_raveler_stitched_grayscale_map(al_l, case_id,  config)
% output_to_raveler_stitched_grayscale(al, config)
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

raveler_config = config.proofreader.Raveler;

save_dir = [get_to_be_proofread_dir(config), raveler_config.grayscale_dir];
check_for_dir(save_dir);

file_name = [save_dir, 'image.', raveler_config.grayscale_version_name, '.', ...
  num2str(case_id, '%05d'), '.png'];
imwrite(uint8(al_l), file_name, 'BitDepth', 8);

return
end
