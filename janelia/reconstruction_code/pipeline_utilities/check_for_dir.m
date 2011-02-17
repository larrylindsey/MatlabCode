function check_for_dir(dir_name)

dir_bracket = strfind(dir_name, '/');
if(isempty(dir_bracket))
  dir_bracket = length(dir_name);
end

for i = 1:length(dir_bracket)
  if(exist(dir_name(1:dir_bracket(i)),'dir')~=7)
    fprintf('making directory %s\n', dir_name);
    mkdir2(dir_name);
    system(sprintf('chmod -R a+rwx "%s"', dir_name));
  end
end

return
end
