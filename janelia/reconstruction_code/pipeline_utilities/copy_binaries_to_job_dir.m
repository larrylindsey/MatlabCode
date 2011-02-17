function copy_binaries_to_job_dir(config, binaries)
% copy_binaries_to_job_dir(binaries)

global config_global

for i = 1:length(binaries)
  copyfile([config_global.bin_dir, binaries{i}], [get_job_dir(config), '.']);
end

return
end
