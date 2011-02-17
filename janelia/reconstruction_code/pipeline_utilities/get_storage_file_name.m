function new_file_name_full = get_storage_file_name(original_file_name_full)
% new_file_name_full = get_storage_file_name(original_file_name_full)
%
% Returns the file name used for storage. If
% config_global.storage.file_name_is_hashed__ is set then file names are
% hashed to reduce length before saving.
%

global config_global

if(~isempty(config_global) && isfield(config_global, 'storage') && ...
    isfield(config_global.storage, 'file_name_is_hashed__') ...
    && config_global.storage.file_name_is_hashed__)
  hash_config = config_global.storage.hash;
  
  LEAVE_AS_IS_PREFIX_LENGTH = hash_config.leave_as_is_prefix_length;
  MIN_LENGTH_TO_HASH = hash_config.min_length_to_hash;
  LEAVE_AS_IS_SUFFIX_LENGTH = hash_config.leave_as_is_suffix_length;
  
  path_length = [0, strfind(original_file_name_full, '/')];
  path_length = path_length(end);

  file_path = original_file_name_full(1:path_length);
  file_name = original_file_name_full(path_length+1:end);

  if(length(file_name)<= (LEAVE_AS_IS_PREFIX_LENGTH + MIN_LENGTH_TO_HASH + ...
      LEAVE_AS_IS_SUFFIX_LENGTH))
    new_file_name = file_name;
  else
    file_name_prefix = file_name(1:LEAVE_AS_IS_PREFIX_LENGTH);
    file_name_suffix = file_name(end-LEAVE_AS_IS_SUFFIX_LENGTH+1:end);
    file_name_to_hash = file_name(LEAVE_AS_IS_PREFIX_LENGTH+1:end-LEAVE_AS_IS_SUFFIX_LENGTH);
    file_name_hashed = hash_string(file_name_to_hash, hash_config);
    new_file_name = [file_name_prefix, file_name_hashed, file_name_suffix];
  end
  new_file_name_full = [file_path, new_file_name];
else
  new_file_name_full = original_file_name_full;
end
return
end
