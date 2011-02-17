function pipeline_cleanup(config)

system(['rm -rf ', config.temp_dir]);

return;
end
