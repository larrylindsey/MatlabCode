function new_file_name_full = save2(original_file_name_full, varargin)
% Save files. If config_file_name_is_hashed__ is set then file names are
% hashed to reduce length before saving.

if(isempty(varargin))
  error('Specify variables to be saved');
end

save_variable_names = ['''', varargin{1}, ''''];
for i = 2:length(varargin)
  save_variable_names = [save_variable_names, ', ''', varargin{i}, ''''];
end

original_file_name = get_file_name_from_full_path(original_file_name_full);

new_file_name_full = get_storage_file_name(original_file_name_full);

execute_exp = '';
if(strcmp(original_file_name_full, new_file_name_full)~=1)
  execute_exp = ['original_file_name_full__ = ''', original_file_name, ''';'];
  save_variable_names = [save_variable_names, ', ''original_file_name_full__'''];
end

load_struct = [];
if(exist(new_file_name_full, 'file')==2)
  warning('off', 'MATLAB:load:variableNotFound');
  load_struct = load(new_file_name_full, 'original_file_name_full__');
  if(isfield(load_struct, 'original_file_name_full__'))
    load_struct.original_file_name_full__ = get_file_name_from_full_path(...
      load_struct.original_file_name_full__);
  end
  warning('on', 'MATLAB:load:variableNotFound');
end
while(~isempty(load_struct) && isfield(load_struct, 'original_file_name_full__') && ...
    strcmp(load_struct.original_file_name_full__, original_file_name)~=1 )
  ext_offset = [length(new_file_name_full)+1, strfind(new_file_name_full, '.')];
  ext_offset = ext_offset(end);
  new_file_name_full = [new_file_name_full(1:ext_offset-1), 'o', ...
    new_file_name_full(ext_offset:end)];
  load_struct = [];
  if(exist(new_file_name_full, 'file')==2)
    warning('off', 'MATLAB:load:variableNotFound');
    load_struct = load(new_file_name_full, 'original_file_name_full__');
    if(isfield(load_struct, 'original_file_name_full__'))
      load_struct.original_file_name_full__ = get_file_name_from_full_path(...
        load_struct.original_file_name_full__);
    end
    warning('on', 'MATLAB:load:variableNotFound');
  end
end

execute_exp = [execute_exp, ...
  '(save(''', new_file_name_full, ''', ', save_variable_names, ', ''-v7.3''))'];
evalin('caller', execute_exp);

return
end
