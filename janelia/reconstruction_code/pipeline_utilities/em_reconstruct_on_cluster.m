function em_reconstruct_on_cluster(config_param, bin_dir, module_ids, ...
  case_ids, job_dir, job_name_prefix, submit_script_name, varargin)
% em_reconstruct_on_cluster(config_param, bin_dir, module_ids, ...
%   case_ids, job_dir, job_name_prefix, submit_script_name, varargin)
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
%   bin_dir             Directory where all the binaries are stored. These
%                         binaries are copied to each job's directory.
%   module_ids          String list of modules to be executed. See pipeline,
%                         pipeline_serial_section, pipeline_block_face for
%                         a list of valid module ids.
%   case_ids            Cell array of sets of sections, will create one job
%                         for each set of sections.
%   job_dir             Directory in which to create the jobs
%   job_name_prefix     String used a prefix for job names, the first
%                         case_id in the job's list is append to this
%                         prefix.
%   submit_script_name  Name of script to store commands for actual
%                         submission to cluster. This script must be run to
%                         perform actual submission to cluster.
%   optional parameter/value list:
%     is_verbose          Whether to print messages [true]
%     is_copied_binary  Whether to copy the binary to each jobs directory
%                         [true].
%     pipeline_bin_name Name of the pipeline binary ['pipeline_sab'].
%     job_file_format   Name format for the jobs. Include %d to print the
%                         first section number for the job's set. E.g, j%d
%                         will produce name j161 for set [161, ...].
%                         Default value is 'j%d'
%     is_emailed_on_error Whether to email the user if the job fails [false]
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  04212009    init. code
%

global config_global code_dir

% Parse the input arguments
is_verbose = true;
is_verbose_figures = false;
is_stand_alone = false;
pipeline_bin_name = 'pipeline_sab';
is_copied_binary = true;
is_emailed_on_error = false;
for i = 1:2:length(varargin)
  switch(varargin{i})
    case 'is_verbose'
      is_verbose = varargin{i+1};
    case 'pipeline_bin_name'
      pipeline_bin_name = varargin{i+1};
    case 'is_copied_binary'
      is_copied_binary = varargin{i+1};
    case 'is_emailed_on_error'
      is_emailed_on_error = varargin{i+1};
    otherwise
      error('Option not understood.');
  end
end

if(is_verbose)
  fprintf('START: em_reconstruct_on_cluster\n');
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
fprintf(fout_submit, ['qsub -N $1 -l excl=true -j y -o log -b y -cwd -V ', ...
  '"%s %s config.mat ''[%s]''"\n'], ...
  wrapper_script, pipeline_bin_name, module_ids);
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
  config.job.name = [job_name_prefix, num2str(case_ids_t(1))];
  config.confirm_email.subject = ['job ', config.job.name, ...
    ': pipeline_sab exited with signal 0'];
  config.confirm_email.text = [config.confirm_email.subject, '\n', ...
    'module: ', module_ids, '\n', ...
    'case_ids: ', num2str(case_ids_t, '%d '), '\n']; 
  
  config_global.is_verbose = is_verbose;
  config_global.is_verbose_figures = is_verbose_figures;
  
  config.region.region_structure = [];
  config.region.region_structure.planes = ...
    region_structure_all.planes(ismember(z_all, case_ids_t));
  
  fprintf('Creating job directories\n');
  job_dir_t = [job_dir, config.job.name, '/'];
  check_for_dir(job_dir_t);
  fprintf('Deleting pre-existing files in job directories\n');
  system(['rm -rf ', job_dir_t, '*']);
  
  fprintf('Saving config file\n');
  save([job_dir_t, 'config.mat'], 'config', 'config_global', 'code_dir');
  
  fprintf('Copying required binaries\n');
  if(is_copied_binary)
    system(['cp -r ', bin_dir, '/* ', job_dir_t, '/.']);
  end
  
  fprintf('Adding commands to submission script\n');
  fprintf(fout_submit, './.%s.submit_job.sh %s\n', job_name_prefix, config.job.name);
end

fclose(fout_submit);
system(['chmod u+x ', job_dir, submit_script_name]);

if(is_verbose)
  fprintf('STOP: em_reconstruct_on_cluster\n');
end

return;
end
