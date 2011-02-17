function seg_eval_prec_recall_bipartite_match(config)
% seg_eval_prec_recall_bipartite_match(config)
% Evaluate segment boundaries by matching with ground-truth
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04282008  init code
% v1 05062008   AG: modified
%

fprintf('Evaluating segment boundaries ..\n');

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
image_dir = get_stack_dir(config);
seg_dir = [get_reconstruction_dir(config), seg_config.dir, evaluate_config.method, '/'];

load2([get_stack_dir(config), evaluate_config.annotations_dir, evaluate_config.groundtruth_file], 'seg');
seg_gt = seg;

evaluation_mode = evaluate_config.mode;

prev_dir = pwd2;
cd(seg_dir);
if(exist(evaluate_config.dir, 'dir')~=7)
  mkdir2(evaluate_config.dir);
end;
cd(prev_dir);

% stack is not prealigned then compute an ROI for the aligned stack if
% needed
if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned ...
    && isempty(stack_config.align.roi))
  align_roi_file_name = [image_dir, stack_config.align.roi_file_name];
  if(exist(align_roi_file_name, 'file')~=2)
    error('Could not find align ROI parameters\n - should have been created during reconstruction\n');
  else
    load2(align_roi_file_name, 'align_roi');
  end;
end

figure_id_prec_recl_plot = 133;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluation of Segment Boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eval_parameters.save_matching = false;
if(strcmp(evaluation_mode, 'match-save')==1)
  eval_parameters.save_matching = true;
end;
eval_parameters.dmax = evaluate_config.max_match_dist; % maximum distance of matching.
eval_parameters.n_images  = size(seg_gt, 3);
eval_parameters.is_verbose = evaluate_config.is_verbose;

precision = [];
recall = [];
match_save_fout = [];

% For dumping a tex for matching results
% Create files for single parameter or a list of parameters
if(isfield(evaluate_config, 'seg_suffix'))
  if(strcmp(evaluation_mode, 'match-save')==1 && isempty(match_save_fout))
    match_save_fout{1} = fopen([seg_dir, evaluate_config.dir, ...
      evaluate_config.master_save_prefix, 'match.', ...
      evaluate_config.seg_suffix, '.tex'], 'w+t');
    create_match_tex(match_save_fout{1}, evaluate_config.seg_suffix);
  end;
else
  if(~iscell(evaluate_config.seg_suffix_id))
    evaluate_config.seg_suffix_id = mat2cell(evaluate_config.seg_suffix_id, ...
      ones(1, size(evaluate_config.seg_suffix_id,1)), ones(1, size(evaluate_config.seg_suffix_id,2)));
  end;
  
  match_save_fout = cell(1, length(evaluate_config.seg_suffix_id));
  for i = 1:length(evaluate_config.seg_suffix_id)
    if(strcmp(evaluation_mode, 'match-save')==1 && (isempty(match_save_fout) || isempty(match_save_fout{i})))
      match_save_fout{i} = fopen(sprintf([seg_dir, evaluate_config.dir, ...
        evaluate_config.master_save_prefix,  'match', ...
        evaluate_config.seg_suffix_format, '.tex'], ...
        evaluate_config.seg_suffix_id{i}), 'w+t');
      create_match_tex(match_save_fout{i}, sprintf(evaluate_config.seg_suffix_format, ...
        evaluate_config.seg_suffix_id{i}));
    end;
  end;
end;

for id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(id);
  fprintf('%d: ', case_id);
  
  image = get_image_from_stack(config, case_id);
  image = image{1};
  if(~isempty(stack_config.roi))
    image = image(stack_config.roi.ymin:stack_config.roi.ymax, ...
      stack_config.roi.xmin:stack_config.roi.xmax);
  end;
  if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned)
    align_config = stack_config.align;
    align_xg_file_name = [image_dir, align_config.xg_file_name];
    xg = read_xf(align_xg_file_name);
    save2([image_dir, align_config.xg_save_name], 'xg');
    
    margin.left = 1-align_config.margin;
    margin.right = size(image,2)+align_config.margin;
    margin.top = 1-align_config.margin;
    margin.bottom = size(image,1)+align_config.margin;
    
    align_tform = xf_2_tform(xg(id,:), size(image, 1), size(image, 2));
    image = apply_tform(image, align_tform, margin, 'nearest');
    image = image(align_roi.ymin:align_roi.ymax, align_roi.xmin:align_roi.xmax);
  end;
  if(evaluate_config.is_verbose_figures)
    figure(1); imshow(image);
  end;
  
  if(isfield(evaluate_config, 'oversegment_ignore_groundtruth_label') && ...
      ~isempty(evaluate_config.oversegment_ignore_groundtruth_label))
    oversegment_ignore_map = ismember(seg_gt(:,:,id), ...
      evaluate_config.oversegment_ignore_groundtruth_label);
  else
    oversegment_ignore_map = [];
  end
  
  segmentation_map_gt = watershed(seg_gt(:,:,id)==0,4);
  if(evaluate_config.is_verbose_figures)
    [py, px] = find(segmentation_map_gt==0);
    figure(1); hold on; plot(px, py, '.'); hold off;
  end
  
  save_prefix = sprintf([seg_dir, evaluate_config.dir, ...
    evaluate_config.master_save_prefix, 'match', '.', ...
    stack_config.image_prefix], case_id);
  
  if(isfield(evaluate_config, 'seg_suffix'))
    segmentation_map = load2(sprintf([seg_dir, stack_config.image_prefix, ...
      evaluate_config.seg_suffix, '.mat'] , case_id),  'label_map');
    
    if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned)
      align_config = stack_config.align;
      align_xg_file_name = [image_dir, align_config.xg_file_name];
      xg = read_xf(align_xg_file_name);
      save2([image_dir, align_config.xg_save_name], 'xg');
      
      margin.left = 1-align_config.margin;
      margin.right = size(segmentation_map.label_map,2)+align_config.margin;
      margin.top = 1-align_config.margin;
      margin.bottom = size(segmentation_map.label_map,1)+align_config.margin;
      
      align_tform = xf_2_tform(xg(id,:), size(segmentation_map.label_map, 1), ...
        size(segmentation_map.label_map, 2));
      segmentation_map.label_map = apply_tform(segmentation_map.label_map, ...
        align_tform, margin, 'nearest');
      segmentation_map.label_map = segmentation_map.label_map(align_roi.ymin:align_roi.ymax, ...
        align_roi.xmin:align_roi.xmax);
    end;
    
    if(isfield(evaluate_config, 'just_display_segmentation') && ...
        ~isempty(evaluate_config.just_display_segmentation) && ...
        evaluate_config.just_display_segmentation)
      figure(2);
      imshow(image);
      [py, px] = find(segmentation_map.label_map==0);
      hold on; plot(px, py, '.'); hold off;
      title(evaluate_config.seg_suffix);
      pause;
      continue;
    end
    
    eval_parameters.save_prefix = [save_prefix, evaluate_config.seg_suffix];
    
    [prec, recl] = evaluate_boundary(segmentation_map.label_map, ...
      segmentation_map_gt, eval_parameters, image, oversegment_ignore_map);
    if(strcmp(evaluation_mode, 'match-save')==1)
      dump_match_tex(match_save_fout{1}, eval_parameters);
    end;
  else
    prec = zeros(length(evaluate_config.seg_suffix_id), 2);
    recl = zeros(length(evaluate_config.seg_suffix_id), 2);
    for i = 1:length(evaluate_config.seg_suffix_id)
      fprintf('%d ', i);
      segmentation_map = load2(sprintf([seg_dir, stack_config.image_prefix, ...
        evaluate_config.seg_suffix_format, '.mat'], case_id, evaluate_config.seg_suffix_id{i}),  'label_map');
      
      if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned)
        align_config = stack_config.align;
        align_xg_file_name = [image_dir, align_config.xg_file_name];
        xg = read_xf(align_xg_file_name);
        save2([image_dir, align_config.xg_save_name], 'xg');
        
        margin.left = 1-align_config.margin;
        margin.right = size(segmentation_map.label_map,2)+align_config.margin;
        margin.top = 1-align_config.margin;
        margin.bottom = size(segmentation_map.label_map,1)+align_config.margin;
        
        align_tform = xf_2_tform(xg(id,:), size(segmentation_map.label_map, 1), ...
          size(segmentation_map.label_map, 2));
        segmentation_map.label_map = apply_tform(segmentation_map.label_map, ...
          align_tform, margin, 'nearest');
        segmentation_map.label_map = segmentation_map.label_map(align_roi.ymin:align_roi.ymax, ...
          align_roi.xmin:align_roi.xmax);
      end;
      
      if(isfield(evaluate_config, 'just_display_segmentation') && ...
          ~isempty(evaluate_config.just_display_segmentation) && ...
          evaluate_config.just_display_segmentation)
        figure(2);
        imshow(image);
        [py, px] = find(segmentation_map.label_map==0);
        hold on; plot(px, py, '.'); hold off;
        title(sprintf(evaluate_config.seg_suffix_format, evaluate_config.seg_suffix_id{i}));
        pause;
        continue;
      end
      
      eval_parameters.save_prefix = sprintf([save_prefix, evaluate_config.seg_suffix_format], ...
        evaluate_config.seg_suffix_id{i});
      
      [prec(i,:), recl(i,:)] = evaluate_boundary(segmentation_map.label_map, ...
        segmentation_map_gt, eval_parameters, image, oversegment_ignore_map);
      if(strcmp(evaluation_mode, 'match-save')==1)
        dump_match_tex(match_save_fout{i}, eval_parameters);
      end;
    end;
  end;
  
  if(isfield(evaluate_config, 'just_display_segmentation') && ...
      ~isempty(evaluate_config.just_display_segmentation) && ...
      evaluate_config.just_display_segmentation)
    continue;
  end
  if(isempty(precision))
    precision = prec;
  else
    precision = precision + prec;
  end;
  if(isempty(recall))
    recall = recl;
  else
    recall = recall + recl;
  end;
  
  
  % Update plot of precision recall curve
  R = recall(:,1) ./ (recall(:,2) + (recall(:,2)==0));
  P = precision(:,1) ./ (precision(:,2) + (precision(:,2)==0));
  figure(figure_id_prec_recl_plot);
  set(gca(), 'FontSize', 16);
  hPlot = plot(R, P, evaluate_config.marker, 'markersize', 12, 'linewidth', 2); hold off;
  title('Segmentation Boundary Evaluation');
  if(isfield(evaluate_config, 'seg_suffix'))
    hLegend = legend(hPlot, evaluate_config.seg_suffix);
    set(hLegend, 'Interpreter', 'none');
  else
    hLegend = legend(hPlot, sprintf(evaluate_config.seg_suffix_format, evaluate_config.seg_suffix_id{end}));
    set(hLegend, 'Interpreter', 'none');
  end;
  axis([0.0 1.0 0.0 1.0]);
  
  fprintf('\n');
end;

if(isfield(evaluate_config, 'just_display_segmentation') && ...
    ~isempty(evaluate_config.just_display_segmentation) && ...
    evaluate_config.just_display_segmentation)
  return;
end

if(strcmp(evaluation_mode, 'pr-curve')==1)
  figure(figure_id_prec_recl_plot);
  if(isfield(evaluate_config, 'seg_suffix'))
    suffix = evaluate_config.seg_suffix;
  else
    suffix = sprintf(evaluate_config.seg_suffix_format, ...
      evaluate_config.seg_suffix_id{end});
  end;
  
  R = recall(:,1) ./ (recall(:,2) + (recall(:,2)==0));
  P = precision(:,1) ./ (precision(:,2) + (precision(:,2)==0));
  [bestR,bestP,bestF] = maxF_2(R,P);
  P_at_R95 = precision_at_given_recall( R, P, .95 );
  P_at_R99 = precision_at_given_recall( R, P, .99 );
  
  write_flag = '';
  if(exist(get_storage_file_name(...
      [seg_dir, evaluate_config.dir, 'seg_boundary_prec_recl', suffix, '.fig']), 'file')==2)
    %     write_flag = input('Overwrite existing plots?([y]/n)', 's');
    write_flag = 'y';
  end
  if(isempty(write_flag) || strcmp(write_flag, 'y')==1)
    saveas(figure_id_prec_recl_plot, ...
      get_storage_file_name([seg_dir, evaluate_config.dir, ...
      evaluate_config.master_save_prefix, 'seg_boundary_prec_recl', ...
      suffix, '.fig']), 'fig');
    save2([seg_dir, evaluate_config.dir, evaluate_config.master_save_prefix, ...
      'seg_boundary_prec_recl', suffix, '.mat'], ...
      'precision', 'recall', 'evaluate_config', 'P_at_R95', 'P_at_R99');
  end;
  
  
  %   system(['chmod a-w ', seg_dir, evaluate_config.dir, 'seg_boundary_*.*']);
end;

if(evaluate_config.is_verbose)
  R = recall(:,1) ./ (recall(:,2) + (recall(:,2)==0));
  P = precision(:,1) ./ (precision(:,2) + (precision(:,2)==0));
  [bestR,bestP,bestF] = maxF_2(R,P);
  P_at_R95 = precision_at_given_recall( R, P, .95 );
  P_at_R99 = precision_at_given_recall( R, P, .99 );
  
  disp(['bestF=' num2str(bestF)]);
  disp(['Precision at 0.95 recall =' num2str(P_at_R95)]);
  disp(['Precision at 0.99 recall =' num2str(P_at_R99)]);
  disp('-- PR eval');
end;

if(strcmp(evaluation_mode, 'match-save')==1)
  for i = 1:length(match_save_fout)
    close_match_tex(match_save_fout{i});
  end;
end;

end

function [prec, recl] = evaluate_boundary(segmentation_map, seg_gt, ...
  eval_parameters, image, oversegment_ignore_map)
% smoothen segmentation map
boundary_map = segmentation_map==0;
boundary_map_d = imdilate(boundary_map, strel('disk', 1));
% boundary_map_d = boundary_map;
segmentation_map = watershed(boundary_map_d,4);

[prec, recl] = evaluate_segmentation_boundary(segmentation_map, seg_gt, ...
  eval_parameters, image, oversegment_ignore_map);
end

function P_at_R = precision_at_given_recall( arrR, arrP, R )
idx=find( arrR>=R, 1,'last' );
if ~isempty(idx)
  if(idx<length(arrR)-1)
    d = (arrR(idx)-R) / (arrR(idx)-arrR(idx+1));
    P_at_R = d*arrP(idx) + (1-d)*arrP(idx+1);
  else
    P_at_R = arrP(idx);
  end
else
  P_at_R =0;
end
end
