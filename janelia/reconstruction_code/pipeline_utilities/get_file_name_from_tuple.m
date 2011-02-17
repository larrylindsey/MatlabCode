function file_name_suffix = get_file_name_from_tuple(save_dir, name_1, name_2, prefix)

file_name_suffix = get_name_prefix_from_string_tuple(name_1, name_2);
dir_bracket = strfind(file_name_suffix, '/');
if(~isempty(dir_bracket))
  dir_bracket = dir_bracket(end);
  check_for_dir([save_dir, file_name_suffix(1:dir_bracket)]);
end
if(nargin>3)
  if(isempty(dir_bracket))
    dir_bracket = 0;
  end
  file_name_suffix = [file_name_suffix(1:dir_bracket), prefix, ...
    file_name_suffix(dir_bracket+1:end)];
end

file_name_suffix = [save_dir, file_name_suffix];

return
end
