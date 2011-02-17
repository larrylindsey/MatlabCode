function name_prefix = get_name_prefix_from_string_tuple(name_1, name_2)
% get_file_prefix_from_string_tuple(name_1, name_2)
% Compute a merged name for two names. Useful for generating file names for
% objects defined over tuples, e.g., transforms between pairs of images.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  01202009  init code
%

global config_global

% get the common part of the path for the names, if any.
dir_brackets_1 = [0, strfind(name_1, '/')];
dir_bracket_1 = dir_brackets_1(end);
subdir_1 = name_1(1:dir_bracket_1);

dir_brackets_2 = [0, strfind(name_2, '/')];
dir_bracket_2 = dir_brackets_2(end);
subdir_2 = name_2(1:dir_bracket_2);

if(~isempty(subdir_1) && ~isempty(subdir_2))
  max_length_common = min(length(subdir_1), length(subdir_2));
  subdir_t_1 = subdir_1(1:max_length_common);
  subdir_t_2 = subdir_2(1:max_length_common);

  is_not_common = abs(subdir_t_1 - subdir_t_2)>0;
  if(nnz(is_not_common)~=0)
    length_common = find(is_not_common, 1) - 1;
  else
    length_common = length(subdir_t_1);
  end
else
  length_common = 0;
end
if(length_common>0)
  length_common = dir_brackets_1(dir_brackets_1<=length_common);
  if(~isempty(length_common))
    length_common = length_common(end);
  else
    length_common = 0;
  end
end

% clean up the names
name_r_1 = strrep(name_1, '/', '_');
name_r_1 = strrep(name_r_1, '\', '_');

name_r_2 = strrep(name_2, '/', '_');
name_r_2 = strrep(name_r_2, '\', '_');

% create a merged name prefix
leave_as_is_prefix_length = ...
  config_global.name_prefix_from_string_tuple.leave_as_is_prefix_length;
name_prefix = [name_1(1:length_common), ...
  name_r_1(dir_bracket_1+1:min(length(name_r_1), dir_bracket_1+leave_as_is_prefix_length)), ...
  '_', ...
  name_r_2(dir_bracket_2+1:min(length(name_r_2), dir_bracket_2+leave_as_is_prefix_length)), ...
  '.', name_r_1, '_', name_r_2];

return
end
