function output_to_raveler_segment_stack(al, seg, output_dir)
% output_to_raveler_segment_stack(al, seg, output_dir)
% Export a segment stack to Raveler so it can be inspected or proofread.
% The grayscale maps are provided in al (Zx1 cell array of MxN images and
% the segment stack is provided in seg (MxNxZ matrix). The files are put in
% the output_dir.
%
% The superpixel maps are the same as segment maps.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

% Dump the grayscale maps
grayscale_dir = [output_dir, 'grayscale_maps/'];
if(exist(grayscale_dir, 'dir')~=7)
  mkdir2(grayscale_dir);
else
  delete([grayscale_dir, '*']);
end
for i = 1:length(al)
  imwrite(al{i}, [grayscale_dir, 'image.v_a', num2str(i, '%05d'), '.png']);
end

% Dump the superpixel maps, superpixel-to-segment-map and
% segment-to-body-map
superpixel_dir = [output_dir, 'superpixel_maps/'];
fout_sp_to_seg = fopen([output_dir, 'superpixel-to-segment-map.txt'], 'wt');
fout_seg_to_body = fopen([output_dir, 'segment-to-body-map.txt'], 'wt');
if(exist(superpixel_dir, 'dir')~=7)
  mkdir2(superpixel_dir);
else
  delete([superpixel_dir, '*']);
end
seg_offset = 0;
for i = 1:size(seg, 3)
  sp = seg(:,:,i);
  imwrite(sp, [superpixel_dir, ...
    'superpixel_map.v_a', num2str(i, '%05d'), '.png'], 'BitDepth', 16);
  
  sp_id = nonzeros(unique(sp(:)));
  seg_id = [0, (1:length(sp_id))+seg_offset];
  sp_id = [0; sp_id]; %#ok<AGROW>
  d = [ones(1, length(sp_id))*i; sp_id', seg_id];
  fprintf(fout_sp_to_seg, '%d %d %d\n', d);
  
  d = [seg_id, sp_id'];
  fprintf(fout_seg_to_body, '%d %d\n', d);
end
fclose(fout_sp_to_seg);
fclose(fout_seg_to_body);

% Dump the 

return
end
