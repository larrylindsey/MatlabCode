function get_full_name(file_name)
try
  l = load(file_name, 'original_file_name_full__');
catch %#ok<CTCH>
  return;
end
if(isempty(l))
  return
end
fprintf('%s\n', l.original_file_name_full__);
return
end
