function mitochondria_segment_2D(config)
% mitochondria_segment_2D(config)
% Segment images to place mitochondria and its mother cell in separate
% segments. Useful for detection.
%
% Nicholas Sofroniew,
% Univ. of Cambridge, UK.
% Visitor March 2009, JFRC, HHMI.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  03112009  init code
%

stack_config = config.stack;
mitochondria_config = config.mitochondria;
segment_config = mitochondria_config.segment;

if(~isfield(mitochondria_config, 'is_verbose'))
  mitochondria_config.is_verbose = true;
end
if(strcmp(segment_config.method, 'marker_controlled_watershed')~=1)
  error('Segment method is not marker controlled watershed');
end

mito_dir = [get_reconstruction_dir(config), mitochondria_config.dir];
check_for_dir(mito_dir);

if(mitochondria_config.is_verbose)
  fprintf('\nSTART: mitochondria_segment_2D\n');
end
for case_id = stack_config.case_ids
  if(mitochondria_config.is_verbose)
    fprintf('plane:%d\n', case_id);
  end
  [images, image_prefixes] = get_image_from_stack(config, case_id);
  image = images{1};
  image_prefix = image_prefixes{1};

  %%%%%%%%%%%
  if(mitochondria_config.is_verbose)
    fprintf('Segmenting ...\n');
  end

  %%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
  %   border=35;
  %   segment_config.minima_depth_threshold=170;
  %   segment_config.hmin_threshold=.01;
  %   segment_config.min_marker_area_threshold=20;
  %   StEl=strel('square',2);
  %   second_minima_size_threshold=50;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if(mitochondria_config.is_verbose)
    fprintf('Computing initial marker mask ... ');
  end
  image = filter_image(image, segment_config.filter_version);
  if(mitochondria_config.is_verbose_figures)
    figure(1);
    imshow(image);
    title('Filtered image');
  end
  image_i = 1 - image;

  % Find image minima with depth below a certain threshold to kill them
  image_i_hmin = imhmin(image_i, segment_config.minima_depth_threshold);
  if(mitochondria_config.is_verbose_figures)
    figure(2);
    imshow(image_i_hmin);
    title('After suppressing minima with depth below a threshold');
  end

  image_i_hmin_i = imsubtract(image_i_hmin, image_i);
  if(mitochondria_config.is_verbose_figures)
    figure(3);
    imshow(image_i_hmin_i);
    title('Inverted hmin image');
  end

  marker_map = im2bw(image_i_hmin_i, segment_config.hmin_threshold);
  if(mitochondria_config.is_verbose_figures)
    figure(4);
    imshow(marker_map);
    title('Markers mask');
  end
  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end

  % Combine minima close to each other so that ideally each mitochondria
  % contains a unique minima
  if(mitochondria_config.is_verbose)
    fprintf('Post processing marker mask ... ');
  end
  marker_map = filter_image(marker_map, segment_config.marker_filter_version);
  if(mitochondria_config.is_verbose_figures)
    figure(5);
    imshow(marker_map);
    title('Markers mask after filtering');
  end
  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end

  % Impose minima on the original image and take watershed
  if(mitochondria_config.is_verbose)
    fprintf('Performing watershed ... ');
  end
  image_i_marked = imimposemin(image_i, marker_map);
  segment_boundary_map_i = logical(watershed(image_i_marked));
  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end

  % Remove mitochondria_config.border and label segments
  if(isfield(mitochondria_config, 'border') && ...
      mitochondria_config.border>0)
    if(mitochondria_config.is_verbose)
      fprintf('Applying border... ');
    end
    segment_boundary_map_i(:,1:mitochondria_config.border-1)=0;
    segment_boundary_map_i(:,end-mitochondria_config.border:end)=0;
    segment_boundary_map_i(1:mitochondria_config.border-1,:)=0;
    segment_boundary_map_i(end-mitochondria_config.border:end,:)=0;
    if(mitochondria_config.is_verbose)
      fprintf('done.\n');
    end
  end

  label_map = bwlabel(segment_boundary_map_i, 8); %#ok<NASGU>
  if(mitochondria_config.is_verbose_figures)
    figure(6);
    imshow(label_map);
    title('Segmentation');
  end
  if(mitochondria_config.is_verbose_figures)
    figure(7);
    imshow(image);
    [py, px] = find(label_map==0);
    hold on;
    plot(px, py, '.');
    hold off;
    title('Segment boundary map');
  end

  if(mitochondria_config.is_verbose)
    fprintf('saving ..');
  end

  save_file_name = [mito_dir, image_prefix, ...
    '.mito.seg.', segment_config.filter_version, ...
    '.', num2str(segment_config.minima_depth_threshold), ...
    '.', num2str(segment_config.hmin_threshold), ...
    '.', segment_config.marker_filter_version, ...
    '.mat'];
  save2(save_file_name, 'label_map');

  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end
end;
if(mitochondria_config.is_verbose)
  fprintf('\nSTOP: mitochondria_segment_2D\n');
end
return;
end
