function output_to_raveler_superpixel_2_seg_map(superpixel_2_seg_map, ...
  superpixel_ids_sets, config)
% output_to_raveler_superpixel_2_seg_map(superpixel_2_seg_map, ...
%   superpixel_ids_sets, config)
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

stack_config = config.stack;
raveler_config = config.proofreader.Raveler;

file_name = [get_to_be_proofread_dir(config), ...
  raveler_config.superpixel_to_segment_file_name];
fout_superpixel_2_seg_map = fopen(file_name, 'wt+');
fprintf(fout_superpixel_2_seg_map, '# sp-segment map, version 0\n');
fprintf(fout_superpixel_2_seg_map, '# plane\tsp\tsegment\n\n');

for plane = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(plane);
  fprintf('plane:%d\n', case_id);
  
  if(~isempty(superpixel_2_seg_map{plane}))
    sp_id = [0, superpixel_ids_sets{plane}]; % 0 is no sp.
    sp_seg_id = superpixel_2_seg_map{plane}(1+sp_id);
    
    d = [repmat(case_id, [1, length(sp_id)]); ...
      sp_id; reshape(sp_seg_id, [1 length(sp_seg_id)])];
    fprintf(fout_superpixel_2_seg_map, '%d\t\t%d\t\t%d\n', d);
  else
    fprintf(fout_superpixel_2_seg_map, '%d\t\t0\t\t0\n', case_id);
  end
  
end
fclose(fout_superpixel_2_seg_map);

return
end
