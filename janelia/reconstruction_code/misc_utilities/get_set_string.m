function set_string = get_set_string(name_strings)
set_string = '.';
for tile_id = 1:length(name_strings)
  set_string = [set_string, name_strings{tile_id}]; %#ok<AGROW>
end
set_string = strrep(set_string, '/', '_');
set_string = strrep(set_string, '\', '_');

return
end
