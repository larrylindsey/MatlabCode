function [annotations_voxel, annotations_body_segment] = get_tile_voxel_annotations(annotations_file_name)
% [annotations_voxel, annotations_body_segment] = get_proofread_voxel_annotations(annotations_file_name)

fin = fopen(annotations_file_name, 'rt');
if(fin==-1)
  error('Could not open annotations file %s\n', annotations_file_name);
end

annotations_voxel = [];
annotations_body_segment = [];
while(~feof(fin))
  annotation_type = fscanf(fin, '%s', 1);
  switch(annotation_type)
    case 'VOXEL'
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

    case 'BODY_SEGMENT'
      a = fgetl(fin);
      % remove characters till the first '{'
      i=1;
      while(strcmp(a(i), '{')~=1 && i<length(a))
        i=i+1;
      end
      % remove characters from the last '}' onwards
      j=length(a);
      while(strcmp(a(j), '}')~=1 && j>1)
        j=j-1;
      end
      if(i==length(a) || j==1 || i+1>=j-1)
        continue;
      end
      segment_ids_list = a(i+1:j-1);
      segment_ids = sscanf(segment_ids_list, '%d');

      % remove characters till the first '"'
      i=1;
      while(strcmp(a(i), '"')~=1 && i<length(a))
        i=i+1;
      end
      a = a(i+1:end);
      while(strcmp(a(end), '"')~=1)
        a = [a, sprintf('\n'), fgetl(fin)]; %#ok<AGROW>
      end
      a = a(1:end-1);

      annob.segments = segment_ids;
      annob.annotation = a;
      
      if(isempty(annotations_body_segment))
        annotations_body_segment = annob;
      else
        annotations_body_segment(end+1) = annob; %#ok<AGROW>
      end
      
    otherwise
      fgetl(fin);
  end
end
fclose(fin);

return
end
