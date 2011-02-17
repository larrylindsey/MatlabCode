function mappings = get_segment_label_remapping(file_name)

fprintf('START: get_segment_label_remapping\n');

mappings = [];
fin = fopen(file_name, 'rt');
while(~feof(fin))
  segment_file_name = fgetl(fin);
  if(isempty(segment_file_name) || strcmp(segment_file_name(1:3), 'Map')~=1)
    continue;
  end
  segment_file_name = segment_file_name(6:end-1);
  fprintf('segment_file_name: %s\n', segment_file_name);
  
  n_map = fgetl(fin);
  if(strcmp(n_map(1:7), 'NUM_MAP')~=1)
    error('keyword NUM_MAP missing');
  end
  n_map = str2double(n_map(9:end));
  fprintf('n_map: %d\n', n_map);
  
  label_mapping = fscanf(fin, '%d', [2 n_map])';
  
  if(isempty(mappings))
    mappings.segment_file_name = segment_file_name;
    mappings.label_mapping = label_mapping;
  else
    mappings(end+1).segment_file_name = segment_file_name; %#ok<AGROW>
    mappings(end).label_mapping = label_mapping; %#ok<AGROW>
  end
end
fprintf('STOP: get_segment_label_remapping\n');
return
end
