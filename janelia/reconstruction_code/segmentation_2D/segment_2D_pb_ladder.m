function segment_2D_pb_ladder(config)
% segment_2D_pb_ladder(config)
% Perform 2D segmentation using pb with ladder
%
% For pb filter see "Learning to Detect Natural Image Boundaries Using
% Local Brightness, Color and Texture Cues" - Martin et al. PAMI04.
%
% Load an image, apply ROI, filter, include mitochondria and vesicles,
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
% v2  04302008  code borrowed from segment_2D_grayscale_ladder.m
% AG: modified 05062008
%

stack_config = config.stack;

seg_config = config.segmentation_2D;
mitochondria_config = config.mitochondria;
vesicle_config = config.vesicle.apply;

if(strcmp(seg_config.method, 'pb_ladder')==0)
  error('Method type does not match with called function. Exiting');
end;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(seg_config.dir, 'dir')~=7)
  mkdir2(seg_config.dir);
end;
cd(seg_config.dir);
if(exist(seg_config.method, 'dir')~=7)
  mkdir2(seg_config.method);
end;
cd(prev_dir);

seg_dir = [get_reconstruction_dir(config), seg_config.dir, seg_config.method, '/'];

pb_model = load2([seg_dir, 'training/', seg_config.model_file, '.mat']);

for case_id = stack_config.case_ids
  fprintf('%d: ', case_id);

  [images, image_prefixes, image_sub_dirs] = get_image_from_stack(config, case_id);

  for tile_id = 1:length(images)
    fprintf('tile %d ', tile_id);
    version_name = 'v0';
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


    %%%%%%%%%%%%%%%%%%%%%%%
    % Compute pb and save for later use
    %%%%%%%%%%%%%%%%%%%%%%%
    version_name = [version_name, '_', seg_config.model_file, '_', ...
      num2str(seg_config.filter_id,'%d_')];
    pb_file = [seg_dir, image_prefix, '.pb.', version_name, '.mat'];
    if(exist(pb_file, 'file')~=2)
      if 1==length(seg_config.filter_id)
        pb = max(applyPBmodelOEOBri( pb_model.model(seg_config.filter_id), ...
          image, seg_config.r_fit_parabola), [],3);
      else
        for is=1:length(seg_config.filter_id)
          pb_s = max(applyPBmodelOEOBri( ...
            pb_model.model(seg_config.filter_id(is)), image), [],3);
          if 1==is
            pb = pb_s;
          else
            pb = max(pb,pb_s);
          end
        end
      end
      min_pb = min(pb(:));
      max_pb = max(pb(:));
      pb = (pb-min_pb)/(max_pb-min_pb);
      save2(pb_file, 'pb');
    else
      load2(pb_file, 'pb');
    end;
    boundary = pb;

    %%%%%%%%%%%%%%%%%%%%%%%
    % Include mitochondria and vesicle if specified
    %%%%%%%%%%%%%%%%%%%%%%%
    if(isfield(seg_config, 'use_mitochondria') && seg_config.use_mitochondria)
      mitochondria_detect = load_mitochondria_detection(...
        [get_reconstruction_dir(config), mitochondria_config.dir, ...
        image_prefix, mitochondria_config.apply.save_suffix, '.mat']);

      mitochondria_mask = mitochondria_detect.mitochondria_detection_confidence > ...
        seg_config.mitochondria_confidence_threshold;

      if(seg_config.mitochondria_erosion>0)
        mitochondria_mask = imerode(mitochondria_mask, strel('disk', seg_config.mitochondria_erosion));
      end;
      if(seg_config.mitochondria_erosion<0)
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


    %%%%%%%%%%%%%%%%%%%%%%%
    % Perform watershed if not already done so
    %%%%%%%%%%%%%%%%%%%%%%%
    check_for_dir([seg_dir, image_sub_dir]);
    watershed_file_name = [seg_dir, image_prefix, '.watershed.', version_name, '.mat'];
    if(exist(watershed_file_name, 'file')==2)
      load2(watershed_file_name,  'watershed_img_c');
    else
      watershed_img_c = watershed(boundary, 4);
      save2(watershed_file_name,  'watershed_img_c'); % direct grayscale
    end;


    %%%%%%%%%%%%%%%%%%%%%%%
    % Compute ladder for the f_thresholds and area_thresholds
    %%%%%%%%%%%%%%%%%%%%%%%
    for s = 1:length(seg_config.f_thresholds)
      for t = 1:length(seg_config.area_thresholds)
        f_threshold = seg_config.f_thresholds(s);
        area_threshold = seg_config.area_thresholds(t);

        save_file_name_suffix = ['.pb_l_T', num2str(f_threshold, '%g'), '_L', num2str(area_threshold, '%g'), '.', version_name];

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
