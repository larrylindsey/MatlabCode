function superpixel_2_segment_grayscale_agglo_mean_boundary(config)
% superpixel_2_segment_grayscale_agglo_mean_boundary(config)
% Perform agglomerative (merger only) segmentation starting from
% superpixels. Criterion: Merge two segments if average boundary value
% between them is less than a threshold
%
% Part of iterative agglomerative clustering methods
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v0  08052008  init code
%

stack_config = config.stack;

seg_config = config.segmentation_2D;
mitochondria_config = config.mitochondria;
vesicle_config = config.vesicle.apply;
superpixel_config = config.superpixel;

if(strcmp(seg_config.method, 'grayscale_AglMeanB')==0)
  error('Method type does not match with called function. Exiting');
end;

seg_dir = get_segmentation_dir(config);
fold_dir = get_fold_dir(config);

for case_id = stack_config.case_ids
  fprintf('%d: ', case_id);
  
  [images, image_prefixes, image_sub_dirs] = get_image_from_stack(config, case_id);
  
  for tile_id = 1:length(images)
    fprintf('tile %d ', tile_id);
    image = images{tile_id};
    image_prefix = image_prefixes{tile_id};
    image_sub_dir = image_sub_dirs{tile_id};
   
    [height, width, n_color] = size(image);
    if(n_color==3)
      image = image(:,:,1);
    end;
    if(~isempty(stack_config.roi))
      image = image(stack_config.roi.ymin:stack_config.roi.ymax, ...
        stack_config.roi.xmin:stack_config.roi.xmax);
    end;
    if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
      image = imresize(image, 1/stack_config.segmentation_scale, ...
      stack_config.segmentation_scaling_method);
    end
    fold_mask = [];
    if(isfield(stack_config, 'fold') && ...
        isfield(stack_config.fold, 'is_considered_in_segmentation') && ...
        stack_config.fold.is_considered_in_segmentation)
      file_name=[fold_dir, image_prefix, '.fold_mask.mat'];
      load2(file_name, 'fold_mask');
      if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
        fold_mask = imresize(fold_mask, 1/stack_config.segmentation_scale, ...
          'nearest');
      end
    end

    [superpixel_method, superpixel_suffixes] = ...
      get_superpixel_suffixes(config, image_prefix);
    superpixel_dir = [get_reconstruction_dir(config), superpixel_config.dir, ...
      superpixel_method, '/'];

    %%%%%%%%%%%%%%%%%%%%%%%
    % Preprocessing stage
    %%%%%%%%%%%%%%%%%%%%%%%
    image = filter_image(image, seg_config.filter_version);
    version_name = seg_config.filter_version;
    if(~isfield(seg_config, 'is_verbose') || seg_config.is_verbose)
      figure(1000); imshow(image);
    end
    
    boundary = 1 - image;


    %%%%%%%%%%%%%%%%%%%%%%%
    % Include mitochondria and vesicle if specified
    %%%%%%%%%%%%%%%%%%%%%%%
    if(~isempty(seg_config.use_mitochondria) && seg_config.use_mitochondria)
      mitochondria_detect = load_mitochondria_detection(...
        [get_reconstruction_dir(config), mitochondria_config.dir, ...
        image_prefix, mitochondria_config.apply.save_suffix, '.mat']);

      mitochondria_mask = mitochondria_detect.mitochondria_detection_confidence>seg_config.mitochondria_confidence_threshold;

      if(seg_config.mitochondria_erosion>0)
        mitochondria_mask = imerode(mitochondria_mask, strel('disk', seg_config.mitochondria_erosion));
      else
        mitochondria_mask = imdilate(mitochondria_mask, strel('disk', -seg_config.mitochondria_erosion));
      end;

      min_area_suffix = '';
      if(isfield(seg_config, 'mitochondria_min_area'))
        mitochondria_mask = bwareaopen(mitochondria_mask, seg_config.mitochondria_min_area);
        min_area_suffix = ['_a', num2str(seg_config.mitochondria_min_area)];
      end;

      if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
        mitochondria_mask = imresize(mitochondria_mask, 1/stack_config.segmentation_scale, ...
          'nearest');
      end
      boundary = boundary .* (1-mitochondria_mask);

      version_name = [version_name, '_m', num2str(seg_config.mitochondria_confidence_threshold), ...
        '_d', num2str(seg_config.mitochondria_erosion), min_area_suffix];
    end;

    if(~isempty(seg_config.use_vesicle) && seg_config.use_vesicle)
      vesicle_detect = load2([get_reconstruction_dir(config), vesicle_config.dir, ...
        image_prefix,vesicle_config.save_suffix, '.mat']);

      vesicle_mask = zeros(size(boundary));
      vesicle_id = find(vesicle_detect.detector_confidence>seg_config.vesicle_threshold);
      vesicle_mask((vesicle_detect.x(vesicle_id)-1)*size(boundary,1) + vesicle_detect.y(vesicle_id))=1;
      vesicle_mask = imdilate(vesicle_mask, strel('disk', seg_config.vesicle_dilation));

      if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
        vesicle_mask = imresize(vesicle_mask, 1/stack_config.segmentation_scale, ...
          'nearest');
      end
      boundary = boundary .* (1-vesicle_mask);

      version_name = [version_name, '_vs'];
    end;

    if(~isempty(fold_mask))
      boundary = max(boundary, fold_mask==0);
      version_name = [version_name, '_f'];
    end

    %%%%%%%%%%%%%%%%%%%%%%%
    % Compute agglomerative clustering for the f_thresholds 
    %%%%%%%%%%%%%%%%%%%%%%%
    fprintf('computing agglomerative clustering: mean boundary .. ');
    check_for_dir([seg_dir, image_sub_dir]);
    
    for superpixel_id = 1:length(superpixel_suffixes)
      superpixel_suffix = superpixel_suffixes{superpixel_id};
      
      superpixel_seg = load2([superpixel_dir, image_prefix, ...
        superpixel_suffix, '.mat']);

      if(~isempty(fold_mask))
        superpixel_seg.label_map = superpixel_seg.label_map .* (fold_mask~=0);
      end

      previous_segment_map = superpixel_seg.label_map;
      if(~isfield(superpixel_seg, 'superpixel_to_seg_label'))
        final_superpixel_2_seg_map = [];
      else
        superpixel_2_seg_map = zeros(max(superpixel_seg.superpixel_to_seg_label(:,1))+1,1);
        superpixel_2_seg_map(superpixel_seg.superpixel_to_seg_label(:,1)+1) = ...
          superpixel_seg.superpixel_to_seg_label(:,2);
        final_superpixel_2_seg_map = superpixel_2_seg_map;
      end

      f_threshold = 0;
      for t = 1:length(seg_config.f_threshold_seq)
        f_threshold_prev = f_threshold;
        f_threshold = seg_config.f_threshold_seq(t);
        length_threshold = seg_config.length_threshold_seq(t);
        fprintf('%g ', f_threshold);
        save_file_name_suffix = ['.gs_amb_sp_T', num2str(f_threshold), ...
          '_', num2str(f_threshold_prev), '_b', num2str(length_threshold), ...
          superpixel_suffix, '.', version_name];

        [label_map, superpixel_to_seg_label] = segmentation_from_superpixels_mean_boundary_value(...
          previous_segment_map, boundary, f_threshold, length_threshold); %#ok<NASGU>

        superpixel_2_seg_map = zeros(max(superpixel_to_seg_label(:,1))+1,1);
        superpixel_2_seg_map(superpixel_to_seg_label(:,1)+1) = ...
          superpixel_to_seg_label(:,2);
        if(isempty(final_superpixel_2_seg_map))
          final_superpixel_2_seg_map = superpixel_2_seg_map;
        else
          final_superpixel_2_seg_map = superpixel_2_seg_map(final_superpixel_2_seg_map+1);
        end
        superpixel_to_seg_label = [find(final_superpixel_2_seg_map)-1, ...
          final_superpixel_2_seg_map(final_superpixel_2_seg_map>0)]; %#ok<NASGU>

        previous_segment_map = label_map;

        segment_boundary = double(label_map==0);
        if(~isfield(seg_config, 'is_verbose') || seg_config.is_verbose)
          figure(2000); imshow(image); hold on;
          [py, px] = find(segment_boundary);
          plot(px, py, '.');
          hold off;
        end;

        if(~isfield(seg_config, 'save_f_thresholds') || ...
            ismember(single(f_threshold), single(seg_config.save_f_thresholds))==1)
          if(~isfield(seg_config, 'superpixel_carry_forward') || ...
              seg_config.superpixel_carry_forward)
            save2([seg_dir, image_prefix, save_file_name_suffix,'.mat'], ...
              'label_map', 'superpixel_to_seg_label');
          else
            save2([seg_dir, image_prefix, save_file_name_suffix,'.mat'], ...
              'label_map');
          end

          save_image = [];
          save_image(:,:,1) = image;
          save_image(:,:,2) = image;
          save_image(:,:,3) = image+segment_boundary;
          imwrite(save_image, ...
            get_storage_file_name([seg_dir, image_prefix, save_file_name_suffix,'.tif']));
        end;
      end;
    end
    fprintf('done\n');
    %   pause;
  end;
end

return
end
