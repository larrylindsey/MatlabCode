function output_stack_annotations(annotations_voxel, ...
  annotations_body, annotations_file_name, config)
% output_stack_annotations(annotations_voxel, ...
%   annotations_body, annotations_file_name, config)

fout = fopen(annotations_file_name, 'wt');
if(fout==-1)
  error('Could not open annotations file %s\n', annotations_file_name);
end

for i = 1:length(annotations_voxel)
  for j = 1:length(annotations_voxel{i})
    fprintf(fout, 'VOXEL %d %d %d "%s"\n', config.stack.case_ids(i), ...
      annotations_voxel{i}(j).y, ...
      annotations_voxel{i}(j).x, annotations_voxel{i}(j).annotation);
  end
end

for i = 1:length(annotations_body)
  for j = 1:length(annotations_body(i).annotation)
    fprintf(fout, 'BODY %d "%s"\n', annotations_body(i).body_id, ...
      annotations_body(i).annotation{j});
  end
end

fclose(fout);

return
end
