function em_reconstruct_on_cluster_import_proofread(config_param, ...
  varargin)
% em_reconstruct_on_cluster_import_proofread(config_param, ...
%   varargin)
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
%   optional
%     is_verbose        Whether to print intermediate messages [true].
%     case_ids          Specify sections to be imported [config.region.case_ids].
%     n_section_per_job  Number of sections per job [1]
%     
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  04282010    init. code
%

% Parse the input arguments
is_verbose = true;
job_dir = '';
job_name_prefix = '';
case_ids = [];
n_section_per_job = 1;
for i = 1:2:length(varargin)
  switch(varargin{i})
    case 'is_verbose'
      is_verbose = varargin{i+1};
    case 'job_dir'
      job_dir = varargin{i+1};
    case 'job_name_prefix'
      job_name_prefix = varargin{i+1};
    case 'case_ids'
      case_ids = varargin{i+1};
    case 'n_section_per_job'
      n_section_per_job = varargin{i+1};
    otherwise
      error('Option not understood.');
  end
end

if(is_verbose)
  fprintf('START: em_reconstruct_on_cluster_import_proofread\n');
end

if(~isempty(job_name_prefix))
  job_name_prefix = [job_name_prefix, '_'];
end

% Get the sections in the region
if(isempty(case_ids))
  if(isa(config_param, 'char'))
    if(strcmp(config_param(end-3:end), '.mat')==1)
      load_config = load(config_param);
      config = load_config.config;
    else
      config = get_basic_config();
      
      h_config_f = str2func(config_param);
      config = h_config_f(config, [], true, false);
    end
  else
    config = config_param.config;
  end
  case_ids = config.region.case_ids;
end

if(isempty(job_dir))
  job_dir = ['cluster_jobs/', num2str(case_ids(1)), '_', ...
    num2str(case_ids(end)), '/'];
end

% Import superpixel maps, adjust them and import superpixel-to-segment maps 
jnp1 = [job_name_prefix, 'is'];
em_reconstruct_on_cluster_synch(config_param, '880.1 880.2 880.3', ...
  get_cell_of_sections(case_ids, 1), job_dir, jnp1);

% Import linkage graphs
jnp2 = [job_name_prefix, 'il'];
em_reconstruct_on_cluster_synch(config_param, '880.4', ...
  get_cell_of_sections(case_ids(1:end-1), 2), job_dir, jnp2, ...
  'prerequisite_job', ['end_', jnp1]);

% Import annotations
jnp3 = [job_name_prefix, 'ia'];
em_reconstruct_on_cluster_synch(config_param, '880.5 880.55', ...
  get_cell_of_sections(case_ids, 1), job_dir, jnp3, ...
  'prerequisite_job', ['end_', jnp1]);

fout = fopen([job_dir, job_name_prefix, 'import.sh'], 'wt');
fprintf(fout, './%sis.sh\n./%sil.sh\n./%sia.sh\n', ...
  job_name_prefix, job_name_prefix, job_name_prefix);
fclose(fout);
system(['chmod u+x ', job_dir, job_name_prefix, 'import.sh']);
fprintf('To submit jobs, run %s%simport.sh\n', job_dir, job_name_prefix);
      
if(is_verbose)
  fprintf('STOP: em_reconstruct_on_cluster_import_proofread\n');
end

return;
end
