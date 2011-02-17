function fout_job = fopen_job_script(config)
% fout_job = fopen_job_script(config)

job_dir = get_job_dir(config);
job_file_name = [job_dir, config.job.name, '.sh'];
if(exist(job_file_name, 'file')~=2)
  error('Job script file has not been initialized correctly. Check prepare_job()');
end

fout_job = fopen(job_file_name, 'at');

return
end
