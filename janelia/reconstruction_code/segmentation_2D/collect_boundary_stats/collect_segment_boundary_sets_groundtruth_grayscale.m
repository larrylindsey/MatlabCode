function collect_segment_boundary_sets_groundtruth_grayscale(config)
% collect_segment_boundary_stats_groundtruth_grayscale(config)
% Collect segment boundary stats from groundtruth.
% Boundary detection algorithm: filtered grayscale.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  08272008  init code
%

train_config = config.train_segmentation;
if(strcmp(train_config.boundary_method, 'grayscale')==0)
  error('Method type does not match with called function. Exiting');
end;

if(~isfield(train_config, 'is_verbose'))
  train_config.is_verbose = true;
end
if(train_config.is_verbose)
  fprintf('Collecting segment boundaries stats from groundtruth ..\n');
end

stack_config = config.stack;
image_dir = get_stack_dir(config);

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
% Collect ground truth segment boundary elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(id);
  if(train_config.is_verbose)
    fprintf('%d ', case_id);
  end
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
  
  segmentation_map_gt = seg_gt(:,:,id); % watershed(seg_gt(:,:,id)==0, 4);
  if(train_config.is_verbose)
    [py, px] = find(segmentation_map_gt==0);
    figure(1); hold on; plot(px, py, '.'); hold off;
  end

  seg_boundary_sets_gt = get_segment_boundary_sets(segmentation_map_gt, boundary_map); %#ok<NASGU>

  check_for_dir([train_dir, image_sub_dir]);
  save([train_dir, image_prefix, '.seg_boundary_sets_gt', '.', ...
    train_config.filter_version, '.mat'], 'seg_boundary_sets_gt');

end
if(train_config.is_verbose)
  fprintf('done\n');
end

return
end
