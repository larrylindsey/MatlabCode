function collect_segment_boundary_sets_auto_segment_grayscale(config)
% collect_segment_boundary_sets_auto_segment(config)
% Collect segment boundary stats from segmentation. The stats are divided
% into correct and over-segmentations.
% Boundary detection algorithm: filtered grayscale.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  08272008  init code
% v1  09022008  code modified for collecting statistics from automatic
%               segmentations.
%

train_config = config.train_segmentation;
if(strcmp(train_config.boundary_method, 'grayscale')==0)
  error('Method type does not match with called function. Exiting');
end;
if(~isfield(train_config, 'is_verbose'))
  train_config.is_verbose = true;
end
if(train_config.is_verbose)
  fprintf('Collecting segment boundaries stats from automatic segmentations ..\n');
end

stack_config = config.stack;
seg_config = config.segmentation_2D;
image_dir = get_stack_dir(config);
seg_dir = [get_reconstruction_dir(config), seg_config.dir, train_config.segment_method, '/'];

load([get_stack_dir(config), train_config.annotations_dir, ...
  train_config.groundtruth_file], 'seg');
seg_gt = seg;

train_dir = get_train_segmentation_dir(config);

% stack is not prealigned then compute an ROI for the aligned stack if
% needed
if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned ...
    && isempty(stack_config.align.roi))
  align_roi_file_name = [image_dir, stack_config.align.roi_file_name];
  if(exist(align_roi_file_name, 'file')~=2)
    error('Could not find align ROI parameters\n - should have been created during reconstruction\n');
  else
    load(align_roi_file_name, 'align_roi');
  end;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collect statistics of Segment Boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eval_parameters.dmax = train_config.max_match_dist; % maximum distance of matching.
eval_parameters.is_verbose = train_config.is_verbose;

if(~isfield(train_config, 'seg_suffix'))
  if(~iscell(train_config.seg_suffix_id))
    train_config.seg_suffix_id = mat2cell(train_config.seg_suffix_id, ...
      ones(1, size(train_config.seg_suffix_id,1)), ones(1, size(train_config.seg_suffix_id,2)));
  end;
end;

for id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(id);
  fprintf('%d: ', case_id);

  [images, image_prefixes, image_sub_dirs] = get_image_from_stack(config, case_id);
  image = images{1};
  image_prefix = image_prefixes{1};
  image_sub_dir = image_sub_dirs{1};
  if(~isempty(stack_config.roi))
    image = image(stack_config.roi.ymin:stack_config.roi.ymax, ...
      stack_config.roi.xmin:stack_config.roi.xmax);
  end;
  if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned)
    align_config = stack_config.align;
    align_xg_file_name = [image_dir, align_config.xg_file_name];
    xg = read_xf(align_xg_file_name);
    save([image_dir, align_config.xg_save_name], 'xg');

    margin.left = 1-align_config.margin;
    margin.right = size(image,2)+align_config.margin;
    margin.top = 1-align_config.margin;
    margin.bottom = size(image,1)+align_config.margin;

    align_tform = xf_2_tform(xg(id,:), size(image, 1), size(image, 2));
    image = apply_tform(image, align_tform, margin, 'nearest');
    image = image(align_roi.ymin:align_roi.ymax, align_roi.xmin:align_roi.xmax);
  end;
  image = filter_image(image, train_config.filter_version);
  if(train_config.is_verbose)
    figure(1); imshow(image);
  end;

  boundary_map = 1 - image;

  if(isfield(train_config, 'oversegment_ignore_groundtruth_label') && ...
      ~isempty(train_config.oversegment_ignore_groundtruth_label))
    oversegment_ignore_map = ismember(seg_gt(:,:,id), ...
      train_config.oversegment_ignore_groundtruth_label);
  else
    oversegment_ignore_map = [];
  end

  segmentation_map_gt = seg_gt(:,:,id); % watershed(seg_gt(:,:,id)==0,4);
  if(train_config.is_verbose)
    [py, px] = find(segmentation_map_gt==0);
    figure(1); hold on; plot(px, py, '.'); hold off;
  end

  check_for_dir([train_dir, image_sub_dir]);
  save_prefix = [train_dir, image_prefix, '.seg_boundary_sets', ...
    '.', train_config.filter_version];

  if(isfield(train_config, 'seg_suffix'))
    segmentation_map = load([seg_dir, image_prefix, train_config.seg_suffix, ...
      '.mat'],  'label_map');

      if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned)
        align_config = stack_config.align;
        align_xg_file_name = [image_dir, align_config.xg_file_name];
        xg = read_xf(align_xg_file_name);
        save([image_dir, align_config.xg_save_name], 'xg');

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
    
    eval_parameters.save_prefix = [save_prefix, train_config.seg_suffix];

    seg_boundary_sets = collect_boundary_sets_auto_segmentation(segmentation_map.label_map, ...
      segmentation_map_gt, eval_parameters, boundary_map, oversegment_ignore_map); %#ok<NASGU>
    save([eval_parameters.save_prefix, '.mat'], 'seg_boundary_sets');
  else
    for i = 1:length(train_config.seg_suffix_id)
      fprintf('%d ', i);
      segmentation_map = load(sprintf([seg_dir, image_prefix, ...
        train_config.seg_suffix_format, '.mat'], train_config.seg_suffix_id{i}), 'label_map');
      
      if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned)
        align_config = stack_config.align;
        align_xg_file_name = [image_dir, align_config.xg_file_name];
        xg = read_xf(align_xg_file_name);
        save([image_dir, align_config.xg_save_name], 'xg');

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

      eval_parameters.save_prefix = sprintf([save_prefix, train_config.seg_suffix_format], ...
        train_config.seg_suffix_id{i});

      seg_boundary_sets = collect_boundary_sets_auto_segmentation(segmentation_map.label_map, ...
        segmentation_map_gt, eval_parameters, boundary_map, oversegment_ignore_map); %#ok<NASGU>
      save([eval_parameters.save_prefix, '.mat'], 'seg_boundary_sets');
    end;
  end;

  if(train_config.is_verbose)
    fprintf('\n');
  end;
end;

if(train_config.is_verbose)
end;

return
end
