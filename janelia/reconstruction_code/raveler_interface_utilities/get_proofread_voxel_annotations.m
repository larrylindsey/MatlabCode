function [annotations_voxel, annotations_body] = get_proofread_voxel_annotations(annotations_file_name)
% annotations_voxel = get_proofread_voxel_annotations(annotations_file_name)

global config_global
if(strcmp(annotations_file_name(end-3:end), '.txt')~=1)
  temp_file_name = [config_global.temp_dir, 'temp_anno.txt'];
  
  system(['$EM_CODE_DIR/raveler_interface_utilities/annotation_pickle_to_ascii.py ', ...
    annotations_file_name, ' > ', temp_file_name]);
  fin = fopen(temp_file_name, 'rt');
  if(fin==-1)
    error('Could not open annotations file %s\n', annotations_file_name);
  end
else
  fin = fopen(annotations_file_name, 'rt');
  if(fin==-1)
    error('Could not open annotations file %s\n', annotations_file_name);
  end
end

annotations_voxel = {};
annotations_body = [];
while(~feof(fin))
  annotation_type = fscanf(fin, '%s', 1);
  switch(annotation_type)
    case 'VOXEL'
      c = fscanf(fin, '%d', [1 3]);
      
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
      
      anno.y = c(2);
      anno.x = c(3);
      anno.annotation = a;
      
      if(c(1)>length(annotations_voxel))
        annotations_voxel{c(1)} = anno; %#ok<AGROW>
      else
        annotations_voxel{c(1)}(end+1) = anno; %#ok<AGROW>
      end
      
    case 'BODY'
      body_id = fscanf(fin, '%d', 1);
      a = fgetl(fin);
      % remove characters till the first '"'
      i=1;
      while(strcmp(a(i), '"')~=1 && i<length(a))
        i=i+1;
      end
      a = a(i+1:end);
      % keep reading lines until a '"' is encountered.
      while(strcmp(a(end),'"')~=1)
        a = [a, sprintf('\n'), fgetl(fin)]; %#ok<AGROW>
      end
      a = a(1:end-1);
      
      if(isempty(annotations_body))
        annotations_body(1).body_id = body_id;
        annotations_body(1).annotations{1} = a;
      else
        a_id = find([annotations_body(:).body_id]==body_id, 1);
        if(~isempty(a_id))
          annotations_body(a_id).annotations{end+1} = a; %#ok<AGROW>
        else
          annotations_body(end+1).body_id = body_id; %#ok<AGROW>
          annotations_body(end).annotations{1} = a; %#ok<AGROW>
        end
      end
    otherwise
      fgetl(fin);
  end
end
fclose(fin);

if(strcmp(annotations_file_name(end-6:end), '.pickle')==1)
  delete(temp_file_name);
end
return
end
