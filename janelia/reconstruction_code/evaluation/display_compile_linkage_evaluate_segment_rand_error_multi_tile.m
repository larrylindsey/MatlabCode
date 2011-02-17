function display_compile_linkage_evaluate_segment_rand_error_multi_tile(config)
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
  fprintf('START: display_compile_linkage_evaluate_segment_rand_error_multi_tile\n');
end

stack_config = config.stack;
linkage_config = config.linkage;
evaluate_config = linkage_config.evaluate;

evaluate_save_dir = [get_region_dir(config), config.linkage.dir, ...
  config.linkage.evaluate.save_dir];

load2([evaluate_save_dir, 'link_eval_segmentwise', ...
  evaluate_config.linkage_suffix_automatic, ...
  evaluate_config.linkage_suffix_ground_truth, '.mat'], 'link_eval');
  
figure;
hPlot = plot(link_eval.rand_score_false_merge, link_eval.rand_score_false_split, ...
  'd-', 'LineWidth', 2, 'MarkerSize', 12);
axis([ 0 max(link_eval.rand_score_false_merge_bound) ...
  0 max(link_eval.rand_score_false_split_bound)]);
title(['Link Eval Segmentwise-Rand ', stack_config.name]);
hLegend = legend(hPlot, [evaluate_config.linkage_suffix_automatic, ...
  evaluate_config.linkage_suffix_ground_truth]);
set(hLegend, 'Interpreter', 'none');
ylabel('Rand Score False Split');
xlabel('Rand Score False Merge');

figure;
hPlot = plot(link_eval.edit_score_n_split, link_eval.edit_score_n_merge, 'd-', ...
  'LineWidth', 2, 'MarkerSize', 12);
axis([0 max(link_eval.n_auto_segment) 0 max(link_eval.n_auto_segment)]);
title(['Link Eval Edit Score ', stack_config.name]);
hLegend = legend(hPlot, [evaluate_config.linkage_suffix_automatic, ...
  evaluate_config.linkage_suffix_ground_truth]);
set(hLegend, 'Interpreter', 'none');
ylabel('Edit Score #Merge');
xlabel('Edit Score #Split');

if(evaluate_config.is_verbose)
  fprintf('STOP: display_compile_linkage_evaluate_segment_rand_error_multi_tile\n');
end
return
end
