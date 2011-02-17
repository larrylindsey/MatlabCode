function compile_boundary_match_score_results(config)
% compile_boundary_match_score_results(config)
% Compile results of evaluation of segment boundaries by matching with
% ground-truth 
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

if(config.is_verbose)
  fprintf('START: compile_boundary_match_score_results\n');
end

evaluate_config = config.eval_segment_boundary;
if(~isfield(evaluate_config,'is_verbose_figures') || ...
    isempty(evaluate_config.is_verbose_figures))
  evaluate_config.is_verbose_figures = true;
end

if(~isfield(evaluate_config, 'master_save_prefix'))
  evaluate_config.master_save_prefix = '';
end

stack_config = config.stack;
seg_config = config.segmentation_2D;
seg_dir = [get_reconstruction_dir(config), seg_config.dir, ...
  evaluate_config.method, '/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compile evaluation results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

precision = zeros(length(evaluate_config.seg_suffix_id),2);
recall = zeros(length(evaluate_config.seg_suffix_id),2);
for id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(id);
  fprintf('case_id: %d\n', case_id);
  
  [image_prefixes, image_sub_dirs] = ...
    get_image_prefixes_subdirs(config, case_id);
  
  for tile_id = 1:length(image_prefixes)
    if(evaluate_config.is_verbose)
      fprintf('tile_id: %d\n', tile_id);
    end
    image_prefix = image_prefixes{tile_id};
    if(evaluate_config.is_verbose)
      fprintf('image_prefix: %s\n', image_prefix);
    end
    image_sub_dir = image_sub_dirs{tile_id};
    save_dir = [seg_dir, evaluate_config.dir, image_sub_dir];
    check_for_dir(save_dir);
    
    load_prefix = [seg_dir, evaluate_config.dir, ...
      evaluate_config.master_save_prefix, 'match', '.', image_prefix];
    
    for i = 1:length(evaluate_config.seg_suffix_id)
      seg_suffix = sprintf(evaluate_config.seg_suffix_format, ...
        evaluate_config.seg_suffix_id{i});
      if(evaluate_config.is_verbose)
        fprintf('seg_suffix: %s\n', seg_suffix);
      end
      load2([load_prefix, seg_suffix, '.', evaluate_config.proofread_name, ...
        '.mat'], 'match_result');
      precision(i,:) = precision(i,:) + match_result.precision;
      recall(i,:) = recall(i,:) + match_result.recall;
    end;
  end;
end

P = precision(:,1) ./ precision(:,2);
R = recall(:,1) ./ recall(:,2);
figure;
set(gca(), 'FontSize', 16);
hPlot = plot(R, P, evaluate_config.marker, 'markersize', 12, 'linewidth', 2);
axis([0 1 0 1]);
title(['Seg Boundary Eval ', stack_config.name, ' ', stack_config.image_prefix]);
hLegend = legend(hPlot, sprintf(evaluate_config.seg_suffix_format, evaluate_config.seg_suffix_id{end}));
set(hLegend, 'Interpreter', 'none');
ylabel('Precision');
xlabel('Recall');

if(config.is_verbose)
  fprintf('STOP: compile_boundary_match_score_results\n');
end
return;
end
