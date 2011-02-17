function em_reconstruct_on_cluster_linkage_train(config_param, train_case_ids, ...
  varargin)
% em_reconstruct_on_cluster_export_for_proofread(config_param)
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
%     job_dir           Directory for the jobs.
%     job_name_prefix   Prefix attached to all job names
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  04282010    init. code
%

% Parse the input arguments
job_dir = ['cluster_jobs/', num2str(train_case_ids(1)), '_', ...
  num2str(train_case_ids(end)), '/'];
job_name_prefix = '';
for i = 1:2:length(varargin)
  switch(varargin{i})
    case 'job_dir'
      job_dir = varargin{i+1};
    case 'job_name_prefix'
      job_name_prefix = varargin{i+1};
    otherwise
      error('Option not understood.');
  end
end

is_verbose = true;
if(is_verbose)
  fprintf('START: em_reconstruct_on_cluster_linkage_train\n');
end

if(~isempty(job_name_prefix))
  job_name_prefix = [job_name_prefix, '_'];
end

% Collect feature vectors
job_name_prefix_full_1 = [job_name_prefix, 'ltf'];
em_reconstruct_on_cluster_synch(config_param, '710.1', ...
  get_cell_of_sections(train_case_ids(1:end-1), 2), job_dir, ...
  job_name_prefix_full_1, 'is_split_by_tiles', true);
    
% Train classifier
job_name_prefix_full_2 = [job_name_prefix, 'ltc'];
em_reconstruct_on_cluster_synch(config_param, '710.5', ...
  {train_case_ids}, job_dir, job_name_prefix_full_2, ...
  'prerequisite_job', ['end_', job_name_prefix_full_1]);

fout = fopen([job_dir, job_name_prefix, 'link_train.sh'], 'wt');
fprintf(fout, './%sltf.sh\n./%sltc.sh\n', ...
  job_name_prefix, job_name_prefix);
fclose(fout);
system(['chmod u+x ', job_dir, job_name_prefix, 'link_train.sh']);
fprintf('To submit jobs, run %s%slink_train.sh\n', job_dir, job_name_prefix);

if(is_verbose)
  fprintf('STOP: em_reconstruct_on_cluster_linkage_train\n');
end

return;
end
