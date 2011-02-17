function job_dir = get_job_dir(config)
% job_dir = get_job_dir(config)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  02282009  init code
%

job_dir = [get_reconstruction_dir(config), config.job.dir, config.job.name, '/'];

return;
end
