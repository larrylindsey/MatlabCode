function mkdir2(new_dir_name)

mkdir(new_dir_name);
if(isunix)
  system(sprintf('chmod -R a+rwx "%s"', new_dir_name));
end

return
end
