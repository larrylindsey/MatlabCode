function output_section_annotations(annotations_voxel_section, ...
  annotations_body_segment_section, annotations_file_name)
% output_section_annotations(annotations_voxel_section, ...
%   annotations_body_segment_section, annotations_file_name)

fout = fopen(annotations_file_name, 'wt');
if(fout==-1)
  error('Could not open annotations file %s\n', annotations_file_name);
end

for i = 1:length(annotations_voxel_section)
%   fprintf('VOXEL %d %d "%s"\n', annotations_voxel_section(i).y, ...
%     annotations_voxel_section(i).x, annotations_voxel_section(i).annotation);
  fprintf(fout, 'VOXEL %d %d "%s"\n', annotations_voxel_section(i).y, ...
    annotations_voxel_section(i).x, annotations_voxel_section(i).annotation);
end

for i = 1:length(annotations_body_segment_section)
%   fprintf('BODY_SEGMENT {');
%   fprintf('%d ', annotations_body_segment_section(i).segments);
%   fprintf('} "%s"\n', annotations_body_segment_section(i).annotation);
  fprintf(fout, 'BODY_SEGMENT {');
  fprintf(fout, '%d ', annotations_body_segment_section(i).segments);
  fprintf(fout, '} "%s"\n', annotations_body_segment_section(i).annotation);
end
fclose(fout);

return
end
