function segment_2D_main(config)
% segment_2D_main(config)
% Perform 2D segmentation. Calls algorithms as specified in config
%
% load an image, apply ROI, filter, include mitochondria and vesicles,
% compute initial watershed, apply ladder and save.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v0  ~Feb.2008 init code
% v1  04112008  modified for reconstruction pipeline
% v2  02112009  modified to unify common code among algorithms
%

stack_config = config.stack;

seg_config = config.segmentation_2D;
mitochondria_config = config.mitochondria;
if(isfield(config, 'vesicle') && isfield(config.vesicle, 'apply'))
  vesicle_config = config.vesicle.apply;
else
  vesicle_config = [];
end

if(~isfield(seg_config, 'is_verbose'))
  seg_config.is_verbose = true;
end
if(~isfield(seg_config, 'is_verbose_figures'))
  seg_config.is_verbose_figures = true;
end
if(seg_config.is_verbose)
  fprintf('START: segment_2D_main\n');
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
    check_for_dir([seg_dir, image_sub_dir]);
    
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
      title('Image after applying ROI and segmentation scale');
    end

    if(isfield(seg_config, 'include_proofread') && ...
        isfield(seg_config.include_proofread, 'is_enabled') && ...
        seg_config.include_proofread.is_enabled)
      proofread_sp_dir = [get_reconstruction_dir(config), ...
        seg_config.dir, config.proofreader.import.dir];
      proofread_sp_suffix = ['.prfrd_asp', '.', ...
        seg_config.include_proofread.proofread_name];
      file_name = [proofread_sp_dir, image_prefix, ...
        proofread_sp_suffix, '.mat'];
      if(seg_config.is_verbose)
        fprintf('Loading proofread sp: %s\n', file_name);
      end
      try
        proofread_sp = load2(file_name, 'label_map', 'locked_labels');
      catch   %#ok<CTCH>
        fprintf('No proofread superpixel map found for this tile.\n');
      end
      proofread_sp.locked_labels = nonzeros(proofread_sp.locked_labels);
      proofread_sp.locked_mask = ismember(proofread_sp.label_map, ...
        proofread_sp.locked_labels);
      proofread_sp.locked_mask = imerode(imdilate(proofread_sp.locked_mask, ...
        strel('disk', 4)), strel('disk', 4));

      if(~isempty(stack_config.roi))
        proofread_sp.locked_mask = proofread_sp.locked_mask(...
          stack_config.roi.ymin:stack_config.roi.ymax, ...
          stack_config.roi.xmin:stack_config.roi.xmax);
        proofread_sp.label_map = proofread_sp.label_map(...
          stack_config.roi.ymin:stack_config.roi.ymax, ...
          stack_config.roi.xmin:stack_config.roi.xmax);
      end;
      if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
        proofread_sp.locked_mask = imresize(proofread_sp.locked_mask, ...
          1/stack_config.segmentation_scale, 'nearest');
        proofread_sp.label_map = imresize(proofread_sp.label_map, ...
          1/stack_config.segmentation_scale, 'nearest');
      end
      if(seg_config.is_verbose)
        fprintf('done.\n');
      end
    else
      proofread_sp = [];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Preprocessing stage
    %%%%%%%%%%%%%%%%%%%%%%%
    if(seg_config.is_verbose)
      fprintf('Computing boundary for watershed ...');
    end
    [boundary_for_watershed, version_name] = compute_boundary_for_watershed();
    if(seg_config.is_verbose)
      fprintf('done\n');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Include mitochondria and vesicle if specified
    %%%%%%%%%%%%%%%%%%%%%%%
    if(seg_config.is_verbose)
      fprintf('Including mitochondria and vesicles if required ...');
    end
    [boundary_for_watershed, version_name] = ...
      include_mitochondria_vesicles_in_boundary(boundary_for_watershed, version_name);
    version_name = regexprep(version_name, '(?<char>[^\\])\', '$<char>');
    if(seg_config.is_verbose)
      fprintf('done\n');
    end
   
    if(seg_config.is_verbose)
      fprintf('Including fold masks if required ...');
    end
    fold_mask = get_fold_mask();
    if(~isempty(fold_mask))
      boundary_for_watershed = max(boundary_for_watershed, fold_mask==0);
      version_name = [version_name, '_f']; %#ok<AGROW>
    end
    if(seg_config.is_verbose)
      fprintf('done\n');
    end

    if(seg_config.is_verbose_figures)
      figure(102);
      imshow(boundary_for_watershed, []);
      title('Boundary map. Intensity rescaled for display.');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Perform watershed if not already done so
    %%%%%%%%%%%%%%%%%%%%%%%
    if(seg_config.is_verbose)
      fprintf('Getting the initial watershed ...');
    end
    if(strcmp(seg_config.method, 'seeded_watershed')==1)
      version_name = regexprep(version_name, '(?<char>[^\\])\', '$<char>');
      watershed_file_name = [seg_dir, image_prefix, '.sd_watershed.', version_name, '.mat'];
      if(exist(watershed_file_name, 'file')~=2)
        watershed_seed_mask = filter_image2(boundary_for_watershed, ...
          seg_config.seeded_watershed.seed_filter, config, image_prefix);
        boundary_for_watershed_seeded = imimposemin(boundary_for_watershed, ...
          watershed_seed_mask, 4);
      else
        boundary_for_watershed_seeded = [];
      end
      watershed_img_c = get_watershed(watershed_file_name, boundary_for_watershed_seeded);
    else
      version_name = regexprep(version_name, '(?<char>[^\\])\', '$<char>');
      watershed_file_name = [seg_dir, image_prefix, '.watershed_img_c.', version_name, '.mat'];
      watershed_img_c = get_watershed(watershed_file_name, boundary_for_watershed);
    end
    if(~isempty(fold_mask))
      watershed_img_c = watershed_img_c .* (fold_mask~=0);
    end
    if(seg_config.is_verbose)
      fprintf('done\n');
    end
    if(isfield(seg_config, 'include_proofread') && ...
        isfield(seg_config.include_proofread, 'is_enabled') && ...
        seg_config.include_proofread.is_enabled)
      if(seg_config.is_verbose)
        fprintf('Removing initial watershed segs in locked regions.\n');
        watershed_img_c(proofread_sp.locked_mask) = 0;
      end
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % compute boundary to be used during segmentation
    %%%%%%%%%%%%%%%%%%%%%%%
    if(seg_config.is_verbose)
      fprintf('Computing boundary for segmentation ...');
    end
    [boundary_for_seg, version_name] = compute_boundary_for_segmentation(version_name);
    version_name = regexprep(version_name, '(?<char>[^\\])\', '$<char>');
    if(seg_config.is_verbose)
      fprintf('done\n');
    end
    if(seg_config.is_verbose_figures)
      figure(103);
      imshow(boundary_for_seg, []);
      title('Boundary map for segmentation. Intensity rescaled for display.');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Call segmentation algo
    %%%%%%%%%%%%%%%%%%%%%%%
    switch(seg_config.method)
      case 'grayscale_ladder'
        segment_2D_sub_ladder(watershed_img_c, boundary_for_seg, image, ...
          seg_config, version_name, [seg_dir, image_prefix], proofread_sp);
      case 'grayscale_AglMeanB'
        segment_2D_sub_agglo_mean_boundary(watershed_img_c, boundary_for_seg, ...
          image, seg_config, version_name, [seg_dir, image_prefix], ...
          proofread_sp);
      case 'seeded_watershed'
        segment_2D_sub_seeded_watershed(watershed_img_c, boundary_for_seg, image, ...
          fold_mask, watershed_seed_mask, seg_config, version_name, [seg_dir, image_prefix], proofread_sp);
      otherwise
        error('Segmentation algorithm type not recognized');
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
  fprintf('STOP: segment_2D_main\n');
end

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

  %%%%%%%%%%%%%%%%%%%%%%%
  % Compute boundary for watershed. A different filter may be used for the
  % segmentation on this watershed.
  %%%%%%%%%%%%%%%%%%%%%%%
  function [boundary_for_watershed, version_name] = compute_boundary_for_watershed()
    if(~isfield(seg_config, 'watershed_filter_version'))
      seg_config.watershed_filter_version = seg_config.filter_version;
    end
    version_name = seg_config.watershed_filter_version;
    image_filtered_watershed = filter_image(...
      image, seg_config.watershed_filter_version, config, image_prefix);
    boundary = 1 - image_filtered_watershed;
    boundary_for_watershed = boundary;
  end

  %%%%%%%%%%%%%%%%%%%%%%%
  % Include mitochondria and vesicle detection into the boundary to reduce
  % over-segmentations.
  %%%%%%%%%%%%%%%%%%%%%%%
  function [b, v] = include_mitochondria_vesicles_in_boundary(b, v)
    if(isfield(seg_config, 'use_mitochondria') && seg_config.use_mitochondria)
      mitochondria_detect = load_mitochondria_detection(...
        [get_reconstruction_dir(config), mitochondria_config.dir, ...
        image_prefix, mitochondria_config.apply.save_suffix, '.mat']);

%       mitochondria_mask = mitochondria_detect.mitochondria_detection_confidence>seg_config.mitochondria_confidence_threshold;
% 
%       if(seg_config.mitochondria_erosion>0)
%         mitochondria_mask = imerode(mitochondria_mask, strel('disk', seg_config.mitochondria_erosion));
%       else
%         mitochondria_mask = imdilate(mitochondria_mask, strel('disk', -seg_config.mitochondria_erosion));
%       end;
% 
%       min_area_suffix = '';
%       if(isfield(seg_config, 'mitochondria_min_area'))
%         mitochondria_mask = bwareaopen(mitochondria_mask, seg_config.mitochondria_min_area);
%         min_area_suffix = ['_a', num2str(seg_config.mitochondria_min_area)];
%       end;
% 
%       if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
%         mitochondria_mask = imresize(mitochondria_mask, 1/stack_config.segmentation_scale, ...
%           'nearest');
%       end
% 
%       % Limit mitochondria mask to areas connected components that aren't
%       % spiculated.  (large perimeter/area metric)
%       mp = regionprops(mitochondria_mask);
%       v = [v, '_m', num2str(seg_config.mitochondria_confidence_threshold), ...
%         '_d', num2str(seg_config.mitochondria_erosion), min_area_suffix];
      
      mt = mitochondria_detect.mitochondria_detection_confidence>seg_config.mitochondria_confidence_threshold;
      mte = imerode(mt, strel('disk', 2));
      mtea = bwareaopen(mte, 500);
      mtead = imdilate(mtea, strel('disk', 2));
      mteadf = imfill(mtead, 'holes');
      mteadfa = bwareaopen(mteadf, 10000);
      l = bwlabel(mteadfa);
      p = regionprops(l, 'Area', 'ConvexArea'); %#ok<MRPBW>
      i = find(([p(:).Area]./[p(:).ConvexArea])<0.85);
      mitochondria_mask = ~ismember(l, [0 i]);
      
      b = b .* (1-mitochondria_mask);

      v = [v, '_m', num2str(seg_config.mitochondria_confidence_threshold), ...
        '_v2'];
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
      b = b .* (1-vesicle_mask);

      v = [v, '_vs'];
    end;
  end

  %%%%%%%%%%%%%%%%%%%%%%%
  % Compute boundary to be used for segmentation.
  %%%%%%%%%%%%%%%%%%%%%%%
  function [boundary, v] = compute_boundary_for_segmentation(v)
    v = [v, '.', seg_config.filter_version];
    image_filtered_seg = filter_image(...
      image, seg_config.filter_version, config, image_prefix);
    boundary = 1 - image_filtered_seg;
    if(~isempty(fold_mask))
      boundary = max(boundary, fold_mask==0);
    end
  end
end

