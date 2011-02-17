function output_to_raveler_stitched_superpixel_map(cat_l, case_id, config)
% output_to_raveler_stitched_superpixel_maps(cat, config)
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

raveler_config = config.proofreader.Raveler;

save_dir = [get_to_be_proofread_dir(config), raveler_config.superpixel_dir];
check_for_dir(save_dir);

if(max(cat_l(:))>2^16)
  error('16 bits is not enough to save the superpixel maps');
end

file_name = [save_dir, 'superpixel_map.', ...
  raveler_config.superpixel_version_name, '.', num2str(case_id, '%05d'), '.png'];
imwrite(uint16(cat_l), file_name, 'BitDepth', 16);

return
end
