function superpixel_2_segment_main(config)
% superpixel_2_segment_main(config)
% Perform 2D segmentation based on superpixels. Calls different
% segmentation algos based on config
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% code borrowed from segment_2D_grayscale_ladder_fly_HRP.m
% v0  01042008  general script
% v1  02212008  fly HRP variant
% v2  03062008  segment from superpixels. Include mitochondria detection,
%               etc.
% v3  04112008  modified to fit into segmentation pipeline framework.
% v4  02112009  modified to unify common code among algorithms
%

stack_config = config.stack;

seg_config = config.segmentation_2D;
mitochondria_config = config.mitochondria;
vesicle_config = config.vesicle.apply;
superpixel_config = config.superpixel;

if(~isfield(seg_config, 'is_verbose'))
  seg_config.is_verbose = true;
end
if(~isfield(seg_config, 'is_verbose_figures'))
  seg_config.is_verbose_figures = true;
end
if(seg_config.is_verbose)
  fprintf('START: superpixel_2_segment_main\n');
end

seg_dir = get_segmentation_dir(config);
fold_dir = get_fold_dir(config);

for case_id = stack_config.case_ids
  if(seg_config.is_verbose)
    fprintf('plane: %d\n', case_id);
  end
  [images, image_prefixes, image_sub_dirs, is_to_be_processed] = ...
    get_image_from_stack(config, case_id);

  for tile_id = 1:length(images)
    image = images{tile_id};
    image_prefix = image_prefixes{tile_id};
    image_sub_dir = image_sub_dirs{tile_id};
    if(seg_config.is_verbose)
      fprintf('tile %d: %s\n', tile_id, image_prefix);
    end
    if(is_to_be_processed(tile_id)==0)
      if(seg_config.is_verbose)
        fprintf('is_to_be_processed=false, skipping\n');
      end
      continue;
    end

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
    if(seg_config.is_verbose_figures)
      figure(101);
      imshow(image);
      title('Image after applying ROI and segmentation_scale');
    end

    [superpixel_method, superpixel_suffixes] = ...
      get_superpixel_suffixes(config, image_prefix);
    superpixel_dir = [get_reconstruction_dir(config), superpixel_config.dir, ...
      superpixel_method, '/'];

    if(isfield(seg_config, 'include_proofread') && ...
        isfield(seg_config.include_proofread, 'is_enabled') && ...
        ~isempty(seg_config.include_proofread.is_enabled) && ...
        seg_config.include_proofread.is_enabled)
      proofread_seg_dir = [get_reconstruction_dir(config), ...
        seg_config.dir, config.proofreader.import.dir];
      proofread_seg_suffix = ['.prfrd_seg', '.', ...
        seg_config.include_proofread.proofread_name];
      file_name = [proofread_seg_dir, image_prefix, ...
        proofread_seg_suffix, '.mat'];
      if(seg_config.is_verbose)
        fprintf('Loading proofread seg: %s\n', file_name);
      end
      proofread_seg = load2(file_name);
      if(seg_config.is_verbose)
        fprintf('done.\n');
      end
    else
      proofread_seg = [];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Preprocessing stage
    %%%%%%%%%%%%%%%%%%%%%%%
    if(seg_config.is_verbose)
      fprintf('Computing boundary for watershed ...');
    end
    [boundary, version_name] = compute_boundary_for_segmentation();
    version_name = regexprep(version_name, '(?<char>[^\\])\', '$<char>');
    if(seg_config.is_verbose)
      fprintf('done\n');
    end


    %%%%%%%%%%%%%%%%%%%%%%%
    % Include mitochondria and vesicle if specified
    %%%%%%%%%%%%%%%%%%%%%%%
    if(seg_config.is_verbose)
      fprintf('Including mitochondria and vesicles if required ...');
    end
    [boundary, version_name, exclusion_mask] = ...
      include_mitochondria_vesicles_in_boundary(boundary, version_name);
    version_name = regexprep(version_name, '(?<char>[^\\])\', '$<char>');
    if(seg_config.is_verbose)
      fprintf('done\n');
    end

    if(seg_config.is_verbose)
      fprintf('Including fold masks if required ...');
    end
    fold_mask = get_fold_mask();
    if(~isempty(fold_mask))
      boundary = max(boundary, fold_mask==0);
      version_name = [version_name, '_f']; %#ok<AGROW>
    end
    if(seg_config.is_verbose)
      fprintf('done\n');
    end

    if(seg_config.is_verbose_figures)
      figure(102);
      imshow(boundary, []);
      title('Boundary map for segmentation. Intensity rescaled for display.');
    end
    
    check_for_dir([seg_dir, image_sub_dir]);
    save_prefix = [seg_dir, image_prefix];
    for suffix_id = 1:length(superpixel_suffixes)
      superpixel_suffix = superpixel_suffixes{suffix_id};

      superpixel_seg = load2([superpixel_dir, image_prefix, ...
        superpixel_suffix, '.mat']);
      if(seg_config.is_verbose)
        fprintf('Superpixel suffix: %s\n', superpixel_suffix);
      end
      superpixel_seg.label_map(superpixel_seg.label_map<0) = 0;
      
      if(~isempty(fold_mask))
        superpixel_seg.label_map = superpixel_seg.label_map .* (fold_mask~=0);
      end

      %%%%%%%%%%%%%%%%%%%%%%%
      % Call segmentation algo
      %%%%%%%%%%%%%%%%%%%%%%%
      switch(seg_config.method)
        case 'grayscale_ladder'
          superpixel_2_segment_sub_ladder(superpixel_seg, boundary, image, ...
            seg_config, version_name, save_prefix, superpixel_suffix, ...
            proofread_seg);
        case {'grayscale_AglBIFL', 'grayscale_AglMeanB', 'grayscale_AglMedianB'}
          superpixel_2_seg_sub_agglo_AMB_ABIFL_AMdB(...
            superpixel_seg, boundary, image, seg_config, version_name, ...
            save_prefix, superpixel_suffix, exclusion_mask, proofread_seg);
        case 'prune_classify'
          superpixel_2_segment_sub_prune_classify(superpixel_seg, ...
            boundary, image, seg_config, version_name, save_prefix, ...
            superpixel_suffix, proofread_seg);
        case 'segment_and'
          superpixel_2_segment_sub_segment_and(superpixel_seg, ...
            boundary, image, seg_config, version_name, save_prefix, ...
            superpixel_suffix, proofread_seg, config, image_prefix);
        otherwise
          error('superpixel-to-segment method not recognized');
      end
    end
    if(seg_config.is_verbose)
      fprintf('done with tile.\n');
    end
  end;
  if(seg_config.is_verbose)
    fprintf('done with plane.\n');
  end
end
if(seg_config.is_verbose)
  fprintf('STOP: superpixel_2_segment_main\n');
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sub routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  % Load the fold mask
  %%%%%%%%%%%%%%%%%%%%%%%
  function fold_mask = get_fold_mask()
    fold_mask = [];
    if(isfield(stack_config, 'fold') && ...
        isfield(stack_config.fold, 'is_considered_in_segmentation') && ...
        stack_config.fold.is_considered_in_segmentation)
      file_name_prefix = [fold_dir, image_prefix, '.fold_mask'];
      fold_mask = load_fold_mask(file_name_prefix, config);
      if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
        fold_mask = imresize(fold_mask, 1/stack_config.segmentation_scale, ...
          'nearest');
      end
    end
  end

  function [boundary, version_name] = compute_boundary_for_segmentation()
    version_name = seg_config.filter_version;
    image_1 = filter_image(...
      image, seg_config.filter_version, config, image_prefix);
    boundary = 1 - image_1;
  end

  function [boundary, version_name, exclusion_mask] = ...
      include_mitochondria_vesicles_in_boundary(boundary, version_name)
    exclusion_mask = zeros(size(image));

    if((isfield(seg_config,'use_mitochondria') && ...
        ~isempty(seg_config.use_mitochondria) && ...
        seg_config.use_mitochondria) || ...
        (isfield(seg_config,'exclude_mitochondria_stat') && ...
        ~isempty(seg_config.exclude_mitochondria_stat) && ...
        seg_config.exclude_mitochondria_stat))
      mitochondria_detect = load_mitochondria_detection(...
        [get_reconstruction_dir(config), mitochondria_config.dir, ...
        image_prefix, mitochondria_config.apply.save_suffix, '.mat']);

      mitochondria_mask = mitochondria_detect.mitochondria_detection_confidence>seg_config.mitochondria_confidence_threshold;

      if(isfield(seg_config,'exclude_mitochondria_stat') && ...
        ~isempty(seg_config.exclude_mitochondria_stat) && ...
          seg_config.exclude_mitochondria_stat)
        exclusion_mask = double(imdilate(mitochondria_mask, strel('disk', 5)));
        if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
          exclusion_mask = imresize(exclusion_mask, 1/stack_config.segmentation_scale, ...
            'nearest');
        end
      end
      
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

      if(isfield(seg_config,'use_mitochondria') && ~isempty(seg_config.use_mitochondria) ...
          && seg_config.use_mitochondria)
        if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
          mitochondria_mask = imresize(mitochondria_mask, 1/stack_config.segmentation_scale, ...
            'nearest');
        end
        boundary = boundary .* (1-mitochondria_mask);

        version_name = [version_name, '_m', num2str(seg_config.mitochondria_confidence_threshold), ...
          '_d', num2str(seg_config.mitochondria_erosion), min_area_suffix];
      end
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
  end

if(seg_config.is_verbose)
  fprintf('STOP: superpixel_2_segment_main\n');
end

return;
end
