function output_matlab_proofread_datastructure_to_raveler(output_dir)

global al cat superpixel_2_seg_map proof

grayscale_dir = [output_dir, 'grayscale_maps/'];
superpixel_map_dir = [output_dir, 'superpixel_maps/'];

% make directories for storage
if(exist(grayscale_dir, 'dir')~=7)
  mkdir(grayscale_dir);
end
if(exist(superpixel_map_dir, 'dir')~=7)
  mkdir(superpixel_map_dir);
end

datastructure_target_dir = output_dir;

% write image files
for plane = 1:length(al)
  file_name = [grayscale_dir, 'image', '.v_a.', num2str(plane, '%05d'), '.png'];
  imwrite(al{plane}, file_name, 'BitDepth', 8);
end

% write superpixel maps
for plane = 1:length(cat)
  file_name = [superpixel_map_dir, ...
    'sp_map', '.v0.', num2str(plane, '%05d'), '.png'];
  imwrite(uint16(cat{plane}), file_name, 'BitDepth', 16);
end

% write superpixel_2_seg_maps
file_name = [datastructure_target_dir, 'superpixel-to-segment-map.txt'];
fout_superpixel_2_seg_map = fopen(file_name, 'wt+');
fprintf(fout_superpixel_2_seg_map, '# sp-segment map\n# plane sp\tsegment\n\n');
for plane = 1:length(superpixel_2_seg_map)
  d = [repmat(plane, [1, length(superpixel_2_seg_map{plane})]); ...
    (1:length(superpixel_2_seg_map{plane}))-1; ...
    reshape(superpixel_2_seg_map{plane}(1:end), [1 length(superpixel_2_seg_map{plane})])];
  fprintf(fout_superpixel_2_seg_map, '%d\t\t%d\t\t%d\n', d);
end
fclose(fout_superpixel_2_seg_map);

% write segment_to_body_map
seg_2_body_map = proof.pmap;

file_name = [datastructure_target_dir, 'segment_to_body_map.txt'];
fout_seg_2_body_map = fopen(file_name, 'wt+');
fprintf(fout_seg_2_body_map, '# segment-body map\n# segment	body\n\n');

d = [(1:length(seg_2_body_map))-1; seg_2_body_map(1:end)];
fprintf(fout_seg_2_body_map, '%d\t\t%d\n', d);

fclose(fout_seg_2_body_map);

return;
end
