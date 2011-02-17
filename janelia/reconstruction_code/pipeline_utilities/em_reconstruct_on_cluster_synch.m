function em_reconstruct_on_cluster_synch(config_param, module_ids, ...
  case_ids, job_dir, job_name_prefix, varargin)
% em_reconstruct_on_cluster_synch(config_param, bin_dir, module_ids, ...
%   case_ids, job_dir, job_name_prefix, varargin)
% executes modules specified in module_id list for the config parameters specified by
% config_param on the cluster
% Inputs
%   config_param        Two options
%                       (1) A '.m' or '.mat' file name. If '.m' then must be a
%                         function returning config. If '.mat', then
%                         variable config is loaded from it.
%                       (2) A struct specifying the reconstruction
%                       parameters through two sub-structs:
%                       config_param.config and config_param.config_global.
%   module_ids          String list of modules to be executed. See pipeline,
%                         pipeline_serial_section, pipeline_block_face for
%                         a list of valid module ids.
%   case_ids            Cell array of sets of sections, will create one job
%                         for each set of sections.
%   job_dir             Directory in which to create the jobs
%   job_name_prefix     String used a prefix for job names, the first
%                         case_id in the job's list is append to this
%                         prefix.
%   optional parameter/value list:
%     is_verbose          Whether to print messages [true]
%     is_copied_binary  Whether to copy the binary to each jobs directory
%                         [true].
%     pipeline_bin_name Name of the pipeline binary ['pipeline_sab'].
%     job_file_format   Name format for the jobs. Include %d to print the
%                         first section number for the job's set. E.g, j_%d
%                         will produce name j161 for set [161, ...].
%                         Default value is 'j_%d'
%     is_emailed_job_done  Whether to email the user when the jobs are
%                         done. This sends a list of failed jobs [true].
%     is_emailed_on_error Whether to email the user if the job fails [false]
%     prerequisite_job  Name of job that must execute before the
%                         to-be-created jobs are allowed to run. Default
%                         empty, i.e. N/A.
%     is_split_by_tiles   Whether to split by tiles in addition to
%                           sections. For each job, a sub-job is created
%                           for the each tile within the first section [false].
%     bin_dir             Directory where all the binaries are stored. These
%                           binaries are copied to each job's directory.
%                           [$EMROOT/bin/pipeline_sab/]
%     submit_script_name  Name of script to store commands for actual
%                           submission to cluster. This script must be run to
%                           perform actual submission to cluster.
%                           [job_name_prefix.sh]
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v1  04292010    stable version
%

global config_global code_dir

% Parse the input arguments
is_verbose = true;
is_verbose_figures = false;
is_stand_alone = false;
pipeline_bin_name = 'pipeline_sab';
is_copied_binary = true;
is_emailed_job_done = true;
is_emailed_on_error = false;
prerequisite_job = '';
is_split_by_tiles = false;
bin_dir = [getenv('EMROOT'), '/bin/pipeline_sab/'];
submit_script_name = [job_name_prefix, '.sh'];
for i = 1:2:length(varargin)
  switch(varargin{i})
    case 'is_verbose'
      is_verbose = varargin{i+1};
    case 'pipeline_bin_name'
      pipeline_bin_name = varargin{i+1};
    case 'is_copied_binary'
      is_copied_binary = varargin{i+1};
    case 'is_emailed_job_done'
      is_emailed_job_done = varargin{i+1};
    case 'is_emailed_on_error'
      is_emailed_on_error = varargin{i+1};
    case 'prerequisite_job'
      prerequisite_job = varargin{i+1};
    case 'is_split_by_tiles'
      is_split_by_tiles = varargin{i+1};
    case 'bin_dir'
      bin_dir = varargin{i+1};
    case 'submit_script_name'
      submit_script_name = varargin{i+1};
    otherwise
      error('Option not understood.');
  end
end

if(is_verbose)
  fprintf('START: em_reconstruct_on_cluster_synch\n');
end

if(~is_copied_binary)
  pipeline_bin_name = [bin_dir, pipeline_bin_name];
else
  pipeline_bin_name = ['./', pipeline_bin_name];
end

wrapper_script = '';
if(is_emailed_on_error)
  wrapper_script = [code_dir, 'lib/misc_utilities/run_email_if_err.sh '];
end

if(~isempty(prerequisite_job))
  prerequisite_job = [' -hold_jid ', prerequisite_job, ' '];
end

if(isa(config_param, 'char'))
  if(strcmp(config_param(end-3:end), '.mat')==1)
    load_config = load(config_param);
    config = load_config.config;
    config_global = load_config.config_global;
  else
    config = get_basic_config();
    
    h_config_f = str2func(config_param);
    config = h_config_f(config, case_ids{1}, is_verbose, is_verbose_figures);
  end
else
  config = config_param.config;
  config_global = config_param.config_global;
end

config_global.job.is_stand_alone = is_stand_alone;
  
config = initialize_region_config_information(config);

% get the region structure: union of all requested case_ids
fprintf('Computing region structure for all sections for all jobs ...\n');
config.stack.case_ids = [];
for i = 1:length(case_ids)
  config.stack.case_ids = [config.stack.case_ids, ...
    case_ids{i}];
end
config.stack.case_ids = unique(config.stack.case_ids);
config_global.temp_dir = '/tmp/';
region_structure_all = get_region_structure(config);
z_all = [region_structure_all.planes(:).z];
fprintf('done.\n');

check_for_dir(job_dir);

fout_submit = fopen([job_dir, '.', job_name_prefix, '.submit_job.sh'], 'wt');
fprintf(fout_submit, 'cd $1\n');
fprintf(fout_submit, ['qsub -N $1 -l excl=true %s -j y -o log -b y -cwd -V ', ...
  '"%s %s config.mat ''[%s]''"\n'], ...
  prerequisite_job, wrapper_script, pipeline_bin_name, module_ids);
fprintf(fout_submit, 'cd ..\n');
fclose(fout_submit);
system(['chmod u+x ', job_dir, '.', job_name_prefix, '.submit_job.sh']);

file_name = [job_dir, submit_script_name];
fprintf('### Job submission script: %s ###\n', file_name);
fout_submit = fopen(file_name, 'wt');
for i = 1:length(case_ids)
  case_ids_t = case_ids{i};
  fprintf('case_ids_t: %d\n', case_ids_t);
  
  fprintf('Constructing config structure\n');
  if(isa(config_param, 'char'))
    if(strcmp(config_param(end-3:end), '.mat')~=1)
      h_config_f = str2func(config_param);
      config = h_config_f(config, case_ids_t, is_verbose, is_verbose_figures);
    else
      config.stack.case_ids = case_ids_t;
    end
  else
    config.stack.case_ids = case_ids_t;
  end
  config_global.job.is_stand_alone = is_stand_alone;
    
  config.is_verbose = is_verbose;
  config.is_verbose_figures = is_verbose_figures;
  config_global.is_verbose = is_verbose;
  config_global.is_verbose_figures = is_verbose_figures;
  
  config.region.region_structure = [];
  config.region.region_structure.planes = ...
    region_structure_all.planes(ismember(z_all, case_ids_t));

  if(~is_split_by_tiles)
    config.job.name = [job_name_prefix, '_', num2str(case_ids_t(1))];
    if(length(config.job.name)>15)
      error('Job name is of length greater than 15. Can cause errors in qsub');
    end
    
    job_dir_t = [job_dir, config.job.name, '/'];
    initialize_job_dir_t(job_dir_t, bin_dir, is_copied_binary);
    
    fprintf('Saving config file\n');
    save([job_dir_t, 'config.mat'], 'config', 'config_global', 'code_dir');
    
    fprintf('Adding commands to submission script\n');
    fprintf(fout_submit, './.%s.submit_job.sh %s\n', job_name_prefix, config.job.name);
  else
    plane_1 = config.region.region_structure.planes(1);
    for j = 1:length(plane_1.tilt_planes(1).tiles)
      plane_1_s = plane_1;
      plane_1_s.tilt_planes.tiles = plane_1.tilt_planes(1).tiles(j);
      config.region.region_structure.planes(1) = plane_1_s;

      config.job.name = [job_name_prefix, '_', num2str(case_ids_t(1)), ...
        '_', num2str(j)];
      if(length(config.job.name)>15)
        error('Job name is of length greater than 15. Can cause errors in qsub');
      end
      
      job_dir_t = [job_dir, config.job.name, '/'];
      initialize_job_dir_t(job_dir_t, bin_dir, is_copied_binary);
      
      fprintf('Saving config file\n');
      save([job_dir_t, 'config.mat'], 'config', 'config_global', 'code_dir');
      
      fprintf('Adding commands to submission script\n');
      fprintf(fout_submit, './.%s.submit_job.sh %s\n', job_name_prefix, config.job.name);
    end
  end
end

fout_end = fopen([job_dir, 'end_', job_name_prefix, '.sh'], 'wt');
fprintf(fout_end, '#!/bin/bash\n');
fprintf(fout_end, 'foo=`ls -1d %s_*/ | chkj | wc | awk ''{print $2;}''`\n', ...
  job_name_prefix);
fprintf(fout_end, 'if [ "$foo" -ne "0" ]\nthen\n');
fprintf(fout_end, 'ls -1d %s_*/ | chkj\n', job_name_prefix);
if(is_emailed_job_done)
  fprintf(fout_end, ['$EM_CODE_DIR/lib/python_email/send_email.sh em.reconstruct ', ...
    'janelia.em.reconstruct $LOGNAME@janelia.hhmi.org %s ', ...
    '"Failed jobs: `ls -1d %s_*/ | chkj | tr [:space:] \\" \\"`"\n'], ...
    job_name_prefix, job_name_prefix);
end
fprintf(fout_end, 'echo "At least one job failed. Deleting all jobs from queue."\nqdel *\nelse\n');
fprintf(fout_end, 'echo "All jobs executed successfully."\n');
if(is_emailed_job_done)
  fprintf(fout_end, ['$EM_CODE_DIR/lib/python_email/send_email.sh em.reconstruct ', ...
    'janelia.em.reconstruct $LOGNAME@janelia.hhmi.org %s "No failed jobs"\nfi\n'], ...
    job_name_prefix);
end
fclose(fout_end);
system(['chmod u+x ', job_dir, 'end_', job_name_prefix, '.sh']);

delete([job_dir, 'end_', job_name_prefix, '_log']);
fprintf(fout_submit, ['qsub -N end_%s -l excl=true -hold_jid "%s_*" -j y ', ...
  '-o end_%s_log -b y -cwd -V ./end_%s.sh\n'], ...
  job_name_prefix, job_name_prefix, job_name_prefix, job_name_prefix);
fclose(fout_submit);
system(['chmod u+x ', job_dir, submit_script_name]);

fprintf('### Job submission script: %s ###\n', [job_dir, submit_script_name]);

if(is_verbose)
  fprintf('STOP: em_reconstruct_on_cluster_synch\n');
end

return;
end

function initialize_job_dir_t(job_dir_t, bin_dir, is_copied_binary)
fprintf('Creating job directories\n');
check_for_dir(job_dir_t);
fprintf('Deleting pre-existing files in job directories\n');
system(['rm -rf ', job_dir_t, '/log']);
fprintf('Copying required binaries\n');
if(is_copied_binary)
  system(['cp -r ', bin_dir, '/* ', job_dir_t, '/.']);
end

return
end
