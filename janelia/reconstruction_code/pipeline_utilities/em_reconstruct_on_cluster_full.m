function em_reconstruct_on_cluster_full(config_param, modes, ...
  varargin)
% em_reconstruct_on_cluster_full(config_param, modes, varargin)
% executes reconstruction modules for the config parameters specified by
% config_param on the cluster
% Inputs
%   config_param        Two options
%                       (1) A '.m' or '.mat' file name. If '.m' then must be a
%                         function returning config. If '.mat', then
%                         variable config is loaded from it.
%                       (2) A struct specifying the reconstruction
%                       parameters through two sub-structs:
%                       config_param.config and config_param.config_global.
%   modes
%     1: superpixel
%     2: superpixel-to-segment
%     3: segment align with section
%     4: 3D linkage
%     5: then run uptil generation of coarse al. These
%     can be viewed to determine proofreading ROI.
%     6: then runs generation of grayscale maps,
%     superpixel maps and other data for proofreading.
%   optional
%     is_verbose        Whether to print intermediate messages [true].
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
for i = 1:2:length(varargin)
  switch(varargin{i})
    case 'is_verbose'
      is_verbose = varargin{i+1};
    case 'job_dir'
      job_dir = varargin{i+1};
    case 'job_name_prefix'
      job_name_prefix = varargin{i+1};
    otherwise
      error('Option not understood.');
  end
end

if(is_verbose)
  fprintf('START: em_reconstruct_on_cluster_export_for_proofread\n');
end

if(~isempty(job_name_prefix))
  job_name_prefix = [job_name_prefix, '_'];
end

% Get the sections in the region
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

if(isempty(job_dir))
  job_dir = ['cluster_jobs/', num2str(case_ids(1)), '_', num2str(case_ids(end)), '/'];
end

for mode = modes
  switch(mode)
    case 1
      % 2D superpixel-to-segment
      jnp1 = [job_name_prefix, 'sp'];
      em_reconstruct_on_cluster_synch(config_param, '450', ...
        get_cell_of_sections(case_ids(1:end-1), 2), job_dir, jnp1);

    case 2
      % 2D superpixel-to-segment
      jnp1 = [job_name_prefix, 'sg'];
      em_reconstruct_on_cluster_synch(config_param, '550', ...
        get_cell_of_sections(case_ids(1:end-1), 2), job_dir, jnp1, ...
        'prerequisite_job', ['end_', job_name_prefix, 'sp']);
      
    case 3
      % For aligning segment maps within sections
      jnp1 = [job_name_prefix, 'as'];
      em_reconstruct_on_cluster_synch(config_param, '650', ...
        get_cell_of_sections(case_ids(1:end-1), 2), job_dir, jnp1, ...
        'prerequisite_job', ['end_', job_name_prefix, 'sg']);
      
    case 4
      % For 3D linkage
      jnp1 = [job_name_prefix, 'l'];
      em_reconstruct_on_cluster_synch(config_param, '750', ...
        get_cell_of_sections(case_ids(1:end-1), 2), job_dir, jnp1, ...
        'prerequisite_job', ['end_', job_name_prefix, 'sg']);
      
    case 5
      % For MATLAB transforms
      jnp1 = [job_name_prefix, 'rt'];
      em_reconstruct_on_cluster_synch(config_param, '810.1', ...
        get_cell_of_sections(case_ids, 1), job_dir, jnp1, ...
        'prerequisite_job', ['end_', job_name_prefix, 'l']);
      
      % For combining MATLAB transforms and for align ROI
      jnp2 = [job_name_prefix, 'rtr'];
      em_reconstruct_on_cluster_synch(config_param, '[810.15, 810.2]', ...
        {case_ids}, job_dir, jnp2, 'prerequisite_job', ['end_', jnp1]);
      
      % For coarse grayscale maps for detemining (manually) the proofreading
      % ROI.
      jnp3 = [job_name_prefix, 'rcg'];
      em_reconstruct_on_cluster_synch(config_param, '810.3', ...
        get_cell_of_sections(case_ids, 1), job_dir, jnp3, ...
        'prerequisite_job', ['end_', jnp2]);
      
      fout = fopen([job_dir, job_name_prefix, 'export_1.sh'], 'wt');
      fprintf(fout, './%srt.sh\n./%srtr.sh\n./%srcg.sh\n', ...
        job_name_prefix, job_name_prefix, job_name_prefix);
      fclose(fout);
      system(['chmod u+x ', job_dir, job_name_prefix, 'export_1.sh']);
      fprintf('To submit jobs, run %s%sexport_1.sh\n', job_dir, job_name_prefix);
      
    case 6
      % For coarse grayscale maps for detemining (manually) the proofreading
      % ROI.
      jnp1 = [job_name_prefix, 'rgs'];
      em_reconstruct_on_cluster_synch(config_param, '[810.35, 810.4]', ...
        get_cell_of_sections(case_ids, 1), job_dir, jnp1, ...
        'prerequisite_job', ['end_', job_name_prefix, 'rtr']);
      
      % For superpixel-to-segment and segment-to-body maps
      jnp2 = [job_name_prefix, 'rsbm'];
      em_reconstruct_on_cluster_synch(config_param, '[810.39, 810.5, 810.6]', ...
        {case_ids}, job_dir, jnp2, 'prerequisite_job', ['end_', jnp1]);
      
      fout = fopen([job_dir, job_name_prefix, 'export_2.sh'], 'wt');
      fprintf(fout, './%srgs.sh\n./%srsbm.sh\n', ...
        job_name_prefix, job_name_prefix);
      fclose(fout);
      system(['chmod u+x ', job_dir, job_name_prefix, 'export_2.sh']);
      fprintf('To submit jobs, run %s%sexport_2.sh\n', job_dir, job_name_prefix);
      
  end
end

if(is_verbose)
  fprintf('STOP: em_reconstruct_on_cluster_export_for_proofread\n');
end

return;
end
