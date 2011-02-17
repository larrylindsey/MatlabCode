function file_name = get_file_name_from_full_path(file_name_full)
% file_name = get_file_name_from_full_path(file_name_full)

dir_bracket = [0, strfind(file_name_full, '/')];
dir_bracket = dir_bracket(end);
file_name = file_name_full(dir_bracket+1:end);

return
end
