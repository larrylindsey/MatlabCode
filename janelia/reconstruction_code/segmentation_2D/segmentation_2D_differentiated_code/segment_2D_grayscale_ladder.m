function segment_2D_grayscale_ladder(config)
% segment_2D_grayscale_ladder(config)
% Perform 2D segmentation using grayscale with ladder
%
% load an image, apply ROI, filter, include mitochondria and vesicles,
% compute initial watershed, apply ladder and save.
%
% Ladder algorithm implemented in
% compute_segmentation_hierarchy_from_watershed_with_min_area_c.cpp
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v0  ~Feb.2008 init code
% v1  04112008  modified for reconstruction pipeline
%

stack_config = config.stack;

seg_config = config.segmentation_2D;
mitochondria_config = config.mitochondria;
vesicle_config = config.vesicle.apply;

if(strcmp(seg_config.method, 'grayscale_ladder')==0)
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

    %%%%%%%%%%%%%%%%%%%%%%%
    % Preprocessing stage
    %%%%%%%%%%%%%%%%%%%%%%%
    version_name = seg_config.filter_version;
    if(isempty(strfind(seg_config.filter_version, ...
        config.normalized_cuts.filter_name_prefix)))
      image = filter_image(image, seg_config.filter_version, config, image_prefix);
%      figure(1000); imshow(image);

      boundary = 1 - image;
    else
      % perform normalized cuts if not already done
      check_for_dir([seg_dir, image_sub_dir]);
      ncut_dir = get_segmentation_ncut_dir(config);
      ncut_file_name = ...
        [ncut_dir, image_prefix, '.', seg_config.filter_version, '.mat'];
      if(exist(ncut_file_name, 'file')==2)
        load2(ncut_file_name,  'boundary');
      else
        boundary = filter_image(image, seg_config.filter_version, config, image_prefix);
        save2(ncut_file_name,  'boundary');
        figure(1000); imshow(boundary, []);
      end;
    end


    %%%%%%%%%%%%%%%%%%%%%%%
    % Include mitochondria and vesicle if specified
    %%%%%%%%%%%%%%%%%%%%%%%
    if(isfield(seg_config, 'use_mitochondria') && seg_config.use_mitochondria)
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
    % Perform watershed if not already done so
    %%%%%%%%%%%%%%%%%%%%%%%
    check_for_dir([seg_dir, image_sub_dir]);
    watershed_file_name = [seg_dir, image_prefix, '.watershed_img_c.', version_name, '.mat'];
    if(exist(watershed_file_name, 'file')==2)
      load2(watershed_file_name,  'watershed_img_c');
    else
      watershed_img_c = watershed(boundary, 4);
      save2(watershed_file_name,  'watershed_img_c'); % direct grayscale
    end;

    if(~isempty(fold_mask))
      watershed_img_c = watershed_img_c .* (fold_mask~=0);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Compute ladder for the f_thresholds and area_thresholds
    %%%%%%%%%%%%%%%%%%%%%%%
    for s = 1:length(seg_config.f_thresholds)
      for t = 1:length(seg_config.area_thresholds)
        f_threshold = seg_config.f_thresholds(s);
        area_threshold = seg_config.area_thresholds(t);

        save_file_name_suffix = ['.gs_l_T', num2str(f_threshold), '_L', num2str(area_threshold), '.', version_name];

        fprintf('computing ladder (hierarchical) .. ');
        label_map = compute_segmentation_hierarchy_from_watershed_with_min_area_c(watershed_img_c, ...
          boundary, f_threshold, area_threshold);
        segment_boundary = double(label_map==0);

        %   fprintf('computing ladder ..');
        %   [r2, r3] = ladder(boundary, thresholds, grayscale_min_area);
        %   segment_boundary = watershed(r3>=thresholds(5), 4)==0;

        figure(2000); imshow(image); hold on;
        [py, px] = find(segment_boundary);
        plot(px, py, '.');
        hold off;

        %     saveas(gcf, sprintf([seg_param.dir, seg_param.image_prefix, save_file_name_suffix,'.fig'] , case_id));
  
        save2([seg_dir, image_prefix, save_file_name_suffix,'.mat'],  'label_map');

        save_image = [];
        save_image(:,:,1) = image;
        save_image(:,:,2) = image;
        save_image(:,:,3) = image+segment_boundary;
        imwrite(save_image, ...
          get_storage_file_name([seg_dir, image_prefix, save_file_name_suffix,'.tif']));
      end;
    end;

    fprintf('done\n');
    %   pause;
  end;
end

end

