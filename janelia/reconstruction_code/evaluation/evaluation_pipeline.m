function evaluation_pipeline(config, main_module_id, sub_module_id)
if(config.is_verbose)
  fprintf('START: evaluation_pipeline\n');
end
switch(main_module_id)
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate segment boundaries using Bipartite matching (CSA++)
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 910
    switch(sub_module_id)
      case 10
        get_boundary_match_score(config);
      case 50
        compile_boundary_match_score_results(config);
    end
  case 911
    error('This is deprecated code. It is advisable to use modules 910.1 and 910.5 for this');
    seg_eval_prec_recall_bipartite_match(config); %#ok<UNRCH>
  case 916
    get_segment_2D_rand_score(config);
  case 917
    get_segment_2D_normalized_rand_score(config);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate linkage confidence estimation using k-fold protocol.
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 920
    switch(sub_module_id)
      case 10
        linkage_evaluate_segment_rand_error(config);
      case 50
        compile_linkage_evaluate_segment_rand_error(config);
      case 51
        display_compile_linkage_evaluate_segment_rand_error(config);
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate reconstruction volume for quantitative analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 930
    prepare_final_reconstruction_volume(config);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate using qscore (qscore3.m)
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 940
    get_q_score(config);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate using Rand (partitionDist.m)
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 941
    get_rand_score(config);

end
if(config.is_verbose)
  fprintf('STOP: evaluation_pipeline\n');
end
return;
end
