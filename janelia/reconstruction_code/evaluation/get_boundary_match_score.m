function get_boundary_match_score(config)
% get_boundary_match_score(config)
% Evaluate segment boundaries by matching with ground-truth
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04282008  init code
% v1 05062008   AG: modified
%

if(config.is_verbose)
  fprintf('START: get_boundary_match_score\n');
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

seg_gt_dir = [get_reconstruction_dir(config), seg_config.dir, ...
  config.proofreader.import.dir];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluation of Segment Boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eval_parameters.dmax = evaluate_config.max_match_dist; % maximum distance of matching.
eval_parameters.is_verbose = evaluate_config.is_verbose;

for id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(id);
  fprintf('case_id: %d\n', case_id);
  
  [images, image_prefixes, image_sub_dirs] = ...
    get_image_from_stack(config, case_id);
  
  for tile_id = 1:length(images)
    if(evaluate_config.is_verbose)
      fprintf('tile_id: %d\n', tile_id);
    end
    image = images{tile_id};
    image_prefix = image_prefixes{tile_id};
    if(evaluate_config.is_verbose)
      fprintf('image_prefix: %s\n', image_prefix);
    end
    image_sub_dir = image_sub_dirs{tile_id};
    save_dir = [seg_dir, evaluate_config.dir, image_sub_dir];
    check_for_dir(save_dir);
    
    if(~isempty(stack_config.roi))
      image = image(stack_config.roi.ymin:stack_config.roi.ymax, ...
        stack_config.roi.xmin:stack_config.roi.xmax);
    end;
    if(evaluate_config.is_verbose_figures)
      figure(1);
      imshow(image);
      title('Image');
    end;
    
    gt = load2([seg_gt_dir, image_prefix, '.prfrd_seg', ...
      '.', evaluate_config.proofread_name, '.mat']);
    seg_gt = gt.label_map;
    if(isfield(evaluate_config, 'oversegment_ignore_groundtruth_label') && ...
        ~isempty(evaluate_config.oversegment_ignore_groundtruth_label))
      oversegment_ignore_map = ismember(seg_gt, ...
        evaluate_config.oversegment_ignore_groundtruth_label);
    else
      oversegment_ignore_map = [];
    end
    
    seg_gt = remove_merged_boundaries_2D(uint32(seg_gt));
    seg_gt = draw_segment_boundaries_2D(uint32(seg_gt));
    segmentation_map_gt = watershed(seg_gt==0,4);
    if(evaluate_config.is_verbose_figures)
      [py, px] = find(segmentation_map_gt==0);
      figure(2);
      imshow(image);
      hold on;
      plot(px, py, '.');
      hold off;
      title('Ground-truth boundaries');
    end
    
    save_prefix = [seg_dir, evaluate_config.dir, ...
      evaluate_config.master_save_prefix, 'match', '.', image_prefix];
    
    for i = 1:length(evaluate_config.seg_suffix_id)
      seg_suffix = sprintf(evaluate_config.seg_suffix_format, ...
        evaluate_config.seg_suffix_id{i});
      if(evaluate_config.is_verbose)
        fprintf('seg_suffix: %s\n', seg_suffix);
      end
      seg = load2([seg_dir, image_prefix, seg_suffix, '.mat'],  'label_map');
      
      % smoothen segmentation map
      boundary_map = seg.label_map==0;
      boundary_map_d = imdilate(boundary_map, strel('disk', 1));
      segmentation_map_as = watershed(boundary_map_d,4);
      
      match_result = evaluate_boundary_bipartite_match(segmentation_map_as==0, ...
        segmentation_map_gt==0, eval_parameters, oversegment_ignore_map); %#ok<NASGU>
      
      save2([save_prefix, seg_suffix, '.', evaluate_config.proofread_name, ...
        '.mat'], 'match_result');
    end;
  end;
end
if(config.is_verbose)
  fprintf('STOP: get_boundary_match_score\n');
end
return;
end
