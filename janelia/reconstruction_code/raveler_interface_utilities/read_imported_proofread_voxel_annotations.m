function annotations_voxel = read_imported_proofread_voxel_annotations(annotations_file_name)
% annotations_voxel = get_proofread_voxel_annotations(annotations_file_name)

fin = fopen(annotations_file_name, 'rt');
if(fin==-1)
  error('Could not open annotations file %s\n', annotations_file_name);
end

annotations_voxel = [];
while(~feof(fin))
  annotation_type = fscanf(fin, '%s', 1);
  if(strcmp(annotation_type, 'VOXEL')==1)
    c = fscanf(fin, '%d', [1 2]);

    a = fgetl(fin);
    % remove characters till the first '"'
    i=1;
    while(strcmp(a(i), '"')~=1 && i<length(a))
      i=i+1;
    end
    % remove characters from the last '"' onwards
    j=length(a);
    while(strcmp(a(j), '"')~=1 && j>1)
      j=j-1;
    end
    if(i==length(a) || j==1 || i+1>=j-1)
      continue;
    end
    a = a(i+1:j-1);
    
    anno.y = c(1);
    anno.x = c(2);
    anno.annotation = a;
    
    if(isempty(annotations_voxel))
      annotations_voxel = anno;
    else
      annotations_voxel(end+1) = anno; %#ok<AGROW>
    end
  else
    fgetl(fin);
  end
end
fclose(fin);
return
end
