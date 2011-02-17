function ex_seg_eval_medulla_on_cluster(mode)
% Example for running segmentation evaluation on the cluster
% Uses a small image from the medulla Leginon dataset.
%
% mode    mode of execution
%         mode = 1 : create job directories.
%         mode = 2 : compile the results into one precision-recall curve
%
% Example of use
% >> ex_seg_eval_medulla_on_cluster(1);
% Log into the cluster and run the job submission scripts.
% >> ex_seg_eval_medulla_on_cluster(2);
%
% Shiv N. Vitaladevuni
%

config = test_seg_eval_medulla(get_basic_config(), 1, true, false);

global config_global

case_ids = 1:2;
seg_suffix_ids = {'0.35_0.34_b0', '0.4_0.39_b0'};

switch(mode)
  case 1
    for i = 1:length(seg_suffix_ids)
      c.config = config;
      c.config.eval_segment_boundary.seg_suffix_id = {seg_suffix_ids{i}};
      c.config_global = config_global;
      em_reconstruct_on_cluster(...
        c, ... % params
        '~/research/em_reconstruction_pipeline/bin/pipeline_sab/', ... % bin dir
        '910.1', ... % module ids
        get_cell_of_sections(case_ids,1), ... % one job per case_id
        'cluster_jobs/test_seg_eval/', ... % job parent dir
        ['se', num2str(i), '_'], ... % job prefix
        ['se', num2str(i), '.sh']); % job submission script name
    end
  case 2
    c.config = config;
    c.config.eval_segment_boundary.seg_suffix_id = seg_suffix_ids;
    c.config_global = config_global;
    em_reconstruct(...
      c, ... % params
      910.5, ... % module ids
      'case_ids', case_ids);
  otherwise
    error('mode not understood');
end

return;
end
