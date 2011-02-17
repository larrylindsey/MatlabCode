function load_result = load2(original_file_name_full, varargin)
% Load files. If config_file_name_is_hashed__ is set then file names are
% hashed to reduce length before loading.

load_variable_names = '';
if(~isempty(varargin))
  load_variable_names = [', ''', varargin{1}, ''''];
  for i = 2:length(varargin)
    load_variable_names = [load_variable_names, ', ''', varargin{i}, ''''];
  end
end

original_file_name = get_file_name_from_full_path(original_file_name_full);

new_file_name_full = get_storage_file_name(original_file_name_full);

is_valid_file = false;
while(exist(new_file_name_full, 'file')==2)
  warning('off', 'MATLAB:load:variableNotFound');
  load_struct = load(new_file_name_full, 'original_file_name_full__');
  warning('on', 'MATLAB:load:variableNotFound');
  if(isfield(load_struct, 'original_file_name_full__'))
    load_struct.original_file_name_full__ = get_file_name_from_full_path(...
      load_struct.original_file_name_full__);
  end
  if(~isfield(load_struct, 'original_file_name_full__') || ...
      strcmp(original_file_name, load_struct.original_file_name_full__)==1)
    is_valid_file = true;
    break;
  end
  ext_offset = [length(new_file_name_full)+1, strfind(new_file_name_full, '.')];
  ext_offset = ext_offset(end);
  new_file_name_full = [new_file_name_full(1:ext_offset-1), 'o', ...
    new_file_name_full(ext_offset:end)];
end
if(~is_valid_file)
  if(exist(original_file_name_full, 'file')==2)
    new_file_name_full = original_file_name_full;
  else
%     fprintf('new_file_name_full: %s\n', new_file_name_full);
%     system(['ls -lt ', new_file_name_full]);
%     fprintf('original_file_name_full: %s\n', original_file_name_full);
%     system(['ls -lt ', original_file_name_full]);
%     dl = [strfind(original_file_name_full, '/'), length(original_file_name_full)];
%     for i = 1:length(dl)
%       fprintf('substr: %s\n', original_file_name_full(1:dl(i)));
%       system(['ls -ld ', original_file_name_full(1:dl(i))]);
%     end
    error(['Could not load file ', original_file_name_full]);
  end
end

if(nargout==0)
  execute_exp = ['(load(''', new_file_name_full, '''', load_variable_names, '))'];
  evalin('caller', execute_exp);
else
  if(isempty(varargin))
    load_result = load(new_file_name_full);
  else
    load_result = load(new_file_name_full, varargin{1});
    for i = 2:length(varargin)
      temp = load(new_file_name_full, varargin{i});
      load_result.(varargin{i})= temp.(varargin{i});
    end
  end
end
return
end
