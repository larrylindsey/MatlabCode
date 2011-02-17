function compile_linkage_evaluate_segment_rand_error_multi_tile(config)
% compile_linkage_evaluate_segment_rand_error_multi_tile(config)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v1  07092008  init code
% v2  06222009  split data section-wise
% v3  07282009  include modified segmentations.
%

if(~isfield(config.linkage.evaluate, 'is_verbose'))
  config.linkage.evaluate.is_verbose = true;
end
if(~isfield(config.linkage.evaluate, 'is_verbose_figures'))
  config.linkage.evaluate.is_verbose_figures = false;
end
if(config.linkage.evaluate.is_verbose)
  fprintf('START: compile_linkage_evaluate_segment_rand_error_multi_tile\n');
end

stack_config = config.stack;
linkage_config = config.linkage;
evaluate_config = linkage_config.evaluate;

evaluate_load_dir = [get_linkage_dir(config), config.linkage.evaluate.save_dir];
evaluate_save_dir = [get_region_dir(config), config.linkage.dir, ...
  config.linkage.evaluate.save_dir];
check_for_dir(evaluate_save_dir);

rand_score_false_split = zeros(length(evaluate_config.linkage_thresholds),1);
rand_score_false_split_bound = zeros(length(evaluate_config.linkage_thresholds),1);
rand_score_false_merge = zeros(length(evaluate_config.linkage_thresholds),1);
rand_score_false_merge_bound = zeros(length(evaluate_config.linkage_thresholds),1);
edit_score_n_split = zeros(length(evaluate_config.linkage_thresholds),1);
edit_score_n_merge = zeros(length(evaluate_config.linkage_thresholds),1);
n_auto_segment = zeros(length(evaluate_config.linkage_thresholds),1);
image_prefixes_2 = {};
for layer_id = 1:length(stack_config.case_ids)-1
  case_id = stack_config.case_ids(layer_id);
  case_id_1 = stack_config.case_ids(layer_id+1);
  fprintf('case_id: %d, case_id_1: %d\n', case_id, case_id_1);
  
  if(isempty(image_prefixes_2))
    image_prefixes_1 = get_image_prefixes_subdirs(config, case_id);
  else
    image_prefixes_1 = image_prefixes_2;
  end
  image_prefixes_2 = get_image_prefixes_subdirs(config, case_id_1);
  
  for tile_1 = 1:length(image_prefixes_1)
    for tile_2 = 1:length(image_prefixes_2)
      if(evaluate_config.is_verbose)
        fprintf('tile 1:%d, tile 2:%d\n', tile_1, tile_2);
        fprintf('tile 1: %s\ntile 2: %s\n', image_prefixes_1{tile_1}, ...
          image_prefixes_2{tile_2});
      end
      
      file_name_prefix = get_file_name_from_tuple(evaluate_load_dir, ...
        image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'elp.');
      file_name_suffix = [evaluate_config.linkage_suffix_automatic, ...
        evaluate_config.linkage_suffix_ground_truth];
      file_name_linkage_eval = [file_name_prefix, file_name_suffix, '.mat'];

      try
        if(evaluate_config.is_verbose)
          fprintf('Attempting to load linkage evaluation:\n%s\n', ...
            file_name_linkage_eval);
        end
        eval_result = load2(file_name_linkage_eval, 'link_eval_p');
      catch %#ok<CTCH>
        if(evaluate_config.is_verbose)
          fprintf('failed, skipping.\n');
        end
        continue;
      end
      if(evaluate_config.is_verbose)
        fprintf('loaded.\n');
      end
      
      linkage_thresholds = ...
        eval_result.link_eval_p.linkage_thresholds;
      rand_score_false_split = rand_score_false_split + ...
        eval_result.link_eval_p.rand_score_false_split;
      rand_score_false_split_bound = rand_score_false_split_bound + ...
        eval_result.link_eval_p.rand_score_false_split_bound;
      rand_score_false_merge = rand_score_false_merge + ...
        eval_result.link_eval_p.rand_score_false_merge;
      rand_score_false_merge_bound = rand_score_false_merge_bound + ...
        eval_result.link_eval_p.rand_score_false_merge_bound;
      edit_score_n_split = edit_score_n_split + ...
        eval_result.link_eval_p.edit_score_n_split;
      edit_score_n_merge = edit_score_n_merge + ...
        eval_result.link_eval_p.edit_score_n_merge;
      n_auto_segment = n_auto_segment + ...
        eval_result.link_eval_p.n_auto_segment;
      
    end
  end
end

link_eval.linkage_thresholds = linkage_thresholds;
link_eval.rand_score_false_split = rand_score_false_split;
link_eval.rand_score_false_split_bound = rand_score_false_split_bound;
link_eval.rand_score_false_merge = rand_score_false_merge;
link_eval.rand_score_false_merge_bound = rand_score_false_merge_bound;
link_eval.edit_score_n_split = edit_score_n_split;
link_eval.edit_score_n_merge = edit_score_n_merge;
link_eval.n_auto_segment = n_auto_segment; %#ok<STRNU>
save2([evaluate_save_dir, 'link_eval_segmentwise', ...
  evaluate_config.linkage_suffix_automatic, ...
  evaluate_config.linkage_suffix_ground_truth, '.mat'], 'link_eval');
  
try
  if(evaluate_config.is_verbose_figures)
    figure;
    hPlot = plot(rand_score_false_merge, rand_score_false_split, 'd-', ...
      'LineWidth', 2, 'MarkerSize', 12);
    axis([ 0 max(rand_score_false_merge_bound) 0 max(rand_score_false_split_bound)]);
    title(['Link Eval Segmentwise-Rand ', stack_config.name]);
    hLegend = legend(hPlot, [evaluate_config.linkage_suffix_automatic, ...
      evaluate_config.linkage_suffix_ground_truth]);
    set(hLegend, 'Interpreter', 'none');
    ylabel('Rand Score False Split');
    xlabel('Rand Score False Merge');
    
    figure;
    hPlot = plot(edit_score_n_split, edit_score_n_merge, 'd-', ...
      'LineWidth', 2, 'MarkerSize', 12);
    axis([0 max(n_auto_segment) 0 max(n_auto_segment)]);
    title(['Link Eval Edit Score ', stack_config.name]);
    hLegend = legend(hPlot, [evaluate_config.linkage_suffix_automatic, ...
      evaluate_config.linkage_suffix_ground_truth]);
    set(hLegend, 'Interpreter', 'none');
    ylabel('Edit Score #Merge');
    xlabel('Edit Score #Split');
  end
catch %#ok<CTCH>
end

if(evaluate_config.is_verbose)
  fprintf('STOP: compile_linkage_evaluate_segment_rand_error_multi_tile\n');
end
return
end
