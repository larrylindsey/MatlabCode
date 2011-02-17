function prepare_job(config)
% prepare_job(config)
%
% Initialize for reconstruction job.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  02292009  init code
%

if(~config.job.is_stand_alone)
  return; % nothing to do
end

job_dir = get_job_dir(config);
check_for_dir(job_dir);

% Preamble
fout = fopen([job_dir, config.job.name, '.sh'], 'wt');
fprintf(fout, '# %s\n', config.job.name);
fprintf(fout, '# Electron Microscope Reconstruction Software\n');
fprintf(fout, '# Chklovskii Lab., Janelia Farm Research Campus,\n');
fprintf(fout, '# Howard Hughes Research Institute.\n');
fprintf(fout, '\n');
%use the 'set -e' option - causes the script to terminate if any run fails
fprintf(fout, 'set -e\n');
% Stack information
fprintf(fout, '# config.stack.name: %s\n', config.stack.name);
if(isfield(config.stack, 'image_structure') && ~isempty(config.stack.image_structure))
  fprintf(fout, '# config.stack.image_structure: %s\n', config.stack.image_structure);
end
fprintf(fout, '# config.stack.case_ids:\n#\t');
for i = 1:length(config.stack.case_ids)
  fprintf(fout, '%d ', config.stack.case_ids(i));
  if(mod(i, 20)==0)
    fprintf(fout, '\n#\t');
  end
end
fprintf(fout, '\n');
fprintf(fout, '# region directory: %s\n', get_region_dir(config));
fprintf(fout, '\n');
% TODO: put config information
fclose(fout);
system(['chmod 777 ', [job_dir, config.job.name, '.sh']]);

return;
end
