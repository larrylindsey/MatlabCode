function get_segment_2D_rand_score(config)
% get_segment_2D_rand_score(config)
% Evaluate segmentation using Rand score by matching with ground-truth
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04282008  init code
% v1  09272009  modified Rand score
%

fprintf('Evaluating segmentation using Rand score ..\n');

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

l = load2([get_stack_dir(config), evaluate_config.annotations_dir, evaluate_config.groundtruth_file], 'seg');
seg_gt = l.seg;

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

if(~iscell(evaluate_config.seg_suffix_id))
  evaluate_config.seg_suffix_id = num2cell(evaluate_config.seg_suffix_id);
end;

figure_id_plot = 133;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluation of Segment Boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eval_parameters.save_matching = false;
eval_parameters.dmax = evaluate_config.max_match_dist; % maximum distance of matching.
eval_parameters.n_images  = size(seg_gt, 3);
eval_parameters.is_verbose = evaluate_config.is_verbose;

rand_score_individual = {};
rand_score_total = [];

for id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(id);
  fprintf('case_id: %d\n', case_id);
  
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
    
    rand_score = evaluate_rand(segmentation_map.label_map, ...
      segmentation_map_gt, oversegment_ignore_map);
  else
    rand_score = zeros(length(evaluate_config.seg_suffix_id), 2);
    for i = 1:length(evaluate_config.seg_suffix_id)
      fprintf('seg. suffix id.: %d\n', i);
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
      
      rand_score(i,:) = evaluate_boundary(segmentation_map.label_map, ...
        segmentation_map_gt, oversegment_ignore_map);
    end;
  end;
  
  if(isfield(evaluate_config, 'just_display_segmentation') && ...
      ~isempty(evaluate_config.just_display_segmentation) && ...
      evaluate_config.just_display_segmentation)
    continue;
  end
  rand_score_individual{id} = rand_score; %#ok<AGROW>
  if(isempty(rand_score_total))
    rand_score_total = rand_score;
  else
    rand_score_total = rand_score_total + rand_score;
  end;
  
%   % Update plot of precision recall curve
%   figure(figure_id_plot);
%   set(gca(), 'FontSize', 16);
%   hPlot = plot(rand_score_total(:,1), rand_score_total(:,2), ...
%     evaluate_config.marker, 'markersize', 12, 'linewidth', 2); hold off;
%   xlabel('rand score false merge');
%   ylabel('rand score false split');
%   title('Segmentation Boundary Evaluation');
%   if(isfield(evaluate_config, 'seg_suffix'))
%     hLegend = legend(hPlot, evaluate_config.seg_suffix);
%     set(hLegend, 'Interpreter', 'none');
%   else
%     hLegend = legend(hPlot, sprintf(evaluate_config.seg_suffix_format, ...
%       evaluate_config.seg_suffix_id{end}));
%     set(hLegend, 'Interpreter', 'none');
%   end;
  
  fprintf('\n');
end;

if(isfield(evaluate_config, 'just_display_segmentation') && ...
    ~isempty(evaluate_config.just_display_segmentation) && ...
    evaluate_config.just_display_segmentation)
  return;
end

% figure(figure_id_plot);
if(isfield(evaluate_config, 'seg_suffix'))
  suffix = evaluate_config.seg_suffix;
else
  suffix = sprintf(evaluate_config.seg_suffix_format, ...
    evaluate_config.seg_suffix_id{end});
end;

write_flag = '';
if(exist(get_storage_file_name(...
    [seg_dir, evaluate_config.dir, 'seg_rand', suffix, '.fig']), 'file')==2)
  %     write_flag = input('Overwrite existing plots?([y]/n)', 's');
  write_flag = 'y';
end
if(isempty(write_flag) || strcmp(write_flag, 'y')==1)
%   saveas(figure_id_plot, ...
%     get_storage_file_name([seg_dir, evaluate_config.dir, ...
%     evaluate_config.master_save_prefix, 'seg_rand', ...
%     suffix, '.fig']), 'fig');
  save2([seg_dir, evaluate_config.dir, evaluate_config.master_save_prefix, ...
    'seg_rand', suffix, '.mat'], 'rand_score_individual', ...
    'rand_score_total', 'evaluate_config');
end;

if(evaluate_config.is_verbose)
  fprintf('best rand score = %g %g\n', ...
    min(rand_score_total(:,1)), min(rand_score_total(:,2)));
end;

end

function rand_score = evaluate_boundary(segmentation_map, seg_gt, ...
  oversegment_ignore_map)
if(~isempty(oversegment_ignore_map))
  seg_gt(oversegment_ignore_map) = 0;
end
seg_gt(seg_gt<0) = 0;
seg_gt = relabel_connected_components_2D(uint32(seg_gt));
rand_score = get_rand_score2(uint32(seg_gt), ...
  uint32(segmentation_map), uint32(10));
end
