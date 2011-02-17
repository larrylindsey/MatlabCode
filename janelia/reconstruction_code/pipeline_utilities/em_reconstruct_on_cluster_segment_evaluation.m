function em_reconstruct_on_cluster_segment_evaluation(config_param, ...
  case_ids, varargin)
% em_reconstruct_on_cluster_segment_evaluation(config_param)
% Inputs
%   config_param        Two options
%                       (1) A '.m' or '.mat' file name. If '.m' then must be a
%                         function returning config. If '.mat', then
%                         variable config is loaded from it.
%                       (2) A struct specifying the reconstruction
%                       parameters through two sub-structs:
%                       config_param.config and config_param.config_global.
%   optional
%     job_dir           Directory for the jobs.
%     job_name_prefix   Prefix attached to all job names
%     segment_modules   Optionally run segmentation modules before
%                         evaluation. Provide list of modules ids.
%     is_split_by_tiles   Whether to split by tiles in addition to
%                           sections. For each job, a sub-job is created
%                           for the each tile within the first section [true].
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  04282010    init. code
%

global config_global

% Parse the input arguments
job_dir = ['cluster_jobs/', num2str(case_ids(1)), '_', ...
  num2str(case_ids(end)), '/'];
job_name_prefix = '';
segment_modules = [];
is_split_by_tiles = true;
for i = 1:2:length(varargin)
  switch(varargin{i})
    case 'job_dir'
      job_dir = varargin{i+1};
    case 'job_name_prefix'
      job_name_prefix = varargin{i+1};
    case 'segment_modules'
      segment_modules = varargin{i+1};
    case 'is_split_by_tiles'
      is_split_by_tiles = varargin{i+1};
    otherwise
      error('Option not understood.');
  end
end

is_verbose = true;
if(is_verbose)
  fprintf('START: em_reconstruct_on_cluster_segment_evaluation\n');
end

if(~isempty(job_name_prefix))
  job_name_prefix = [job_name_prefix, '_'];
end

if(isa(config_param, 'char'))
  if(strcmp(config_param(end-3:end), '.mat')==1)
    load_config = load(config_param);
    config = load_config.config;
    config_global = load_config.config_global;
  else
    config = get_basic_config();
    
    h_config_f = str2func(config_param);
    config = h_config_f(config, case_ids, true, false);
  end
else
  config = config_param.config;
  config_global = config_param.config_global;
end

fout = fopen([job_dir, job_name_prefix, 'eval_seg.sh'], 'wt');

if(~isempty(segment_modules))
  job_name_prefix_full = [job_name_prefix, 's'];
  em_reconstruct_on_cluster_synch(config_param, num2str(segment_modules, '%g '), ...
    get_cell_of_sections(case_ids, 1), job_dir, ...
    job_name_prefix_full, 'is_split_by_tiles', is_split_by_tiles);
  fprintf(fout, './%s.sh\n', job_name_prefix_full);
end

seg_eval_suffixes = config.eval_segment_boundary.seg_suffix_id;
for t = 1:length(seg_eval_suffixes)
  config_param_t.config = config;
  config_param_t.config_global = config_global;
  config_param_t.config.eval_segment_boundary.seg_suffix_id = ...
    config.eval_segment_boundary.seg_suffix_id(t);
  
  job_name_prefix_full_1 = [job_name_prefix, 'es_', num2str(t)];
  if(isempty(segment_modules))
    em_reconstruct_on_cluster_synch(config_param_t, '910.1', ...
      get_cell_of_sections(case_ids, 1), job_dir, ...
      job_name_prefix_full_1, 'is_split_by_tiles', is_split_by_tiles);
  else
    em_reconstruct_on_cluster_synch(config_param_t, '910.1', ...
      get_cell_of_sections(case_ids, 1), job_dir, ...
      job_name_prefix_full_1, 'is_split_by_tiles', is_split_by_tiles, ...
      'prerequisite_job', [job_name_prefix, 'end_s']);
  end
  fprintf(fout, './%s.sh\n', job_name_prefix_full_1);
end
fclose(fout);
system(['chmod u+x ', job_dir, job_name_prefix, 'eval_seg.sh']);
fprintf('To submit jobs, run %s%seval_seg.sh\n', job_dir, job_name_prefix);

if(is_verbose)
  fprintf('STOP: em_reconstruct_on_cluster_segment_evaluation\n');
end

return;
end
