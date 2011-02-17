function generate_annotations_multi_tile_patches(config)
% generate_annotations_multi_tile_patches(config)
% generate annotations when there are multiple tiles and patches per section (
% trakEM xml)
%
% al{z}'s are the EM images used for the reconstruction
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
% v1  09302008  modified to multiple tile per section (trakEM xml)
%

if(config.is_verbose)
  fprintf('START: generate_annotations_multi_tile_patches\n');
end
stack_config = config.stack;

switch(config.align.global_align.method)
  case 'SIFT'
    tform_dir = [get_region_dir(config), config.SIFT.dir];
  case 'deformable_mesh'
    tform_dir = [get_region_dir(config), config.deformable_mesh.dir];
end

load2([get_global_stitching_param_dir(config), 'canvas.mat'], ...
  'canvas_height', 'canvas_width');

align_roi_file_name = [get_global_stitching_param_dir(config), ...
  config.stack.align.roi_file_name];
load2(align_roi_file_name, 'align_roi');

annot_dir = [get_reconstruction_dir(config), ...
  config.annotations.dir, config.proofreader.annotations.method, '/'];
fold_dir = get_fold_dir(config);
align_seg_dir = get_align_seg_dir(config);

for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  fprintf('plane: %d\n', case_id);
  [image_prefixes, image_sub_dirs] = ...
    get_image_prefixes_subdirs(config, case_id);
  
  load2([get_global_stitching_param_dir(config), ...
    'global_stitching_parameters.c.', ...
    num2str(case_id), '.mat'], 'xdata', 'ydata', 'tforms');
  
  % load label mapping
  image_set_string = get_set_string(image_prefixes);
  align_seg_save_dir = [align_seg_dir, image_sub_dirs{1}];
  label_mapping = load2([align_seg_save_dir, 'sec.', num2str(case_id), ...
    '.aligned_seg_mapping', image_set_string, ...
    config.segmentation_choose.choice.seg_suffix, '.mat'], 'label_mappings');
  for tile = 1:length(image_prefixes)
    label_mapping.label_mappings{tile}(label_mapping.label_mappings{tile}<0)=0;
  end
  
  annotations_voxel_section = [];
  annotations_body_segment_section = [];
  annotations_file_name = [get_region_dir(config), config.annotations.dir, ...
    'annotations.', num2str(case_id), '.txt'];
  delete(annotations_file_name);
  
  section_image = zeros(canvas_height, canvas_width)-1;
  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))
      continue;
    end
    fprintf('tile: %d, %s\n', tile_id, image_prefixes{tile_id});
    
    transform_affine = {};
    if(~isfield(config.align.global_align, 'storage_format') || ...
        strcmp(config.align.global_align.storage_format, 'mat')==1)
      load2([tform_dir, image_prefixes{tile_id}, ...
        '.patchwise_affine_transforms.mat'], 'transform');
      n_patch = size(transform.transform, 2);
      for patch_id = 1:n_patch
        if(config.is_verbose)
          fprintf('patch: %d\n', patch_id);
        end
        transform_affine{patch_id} = ...
          reshape(transform.transform(:,patch_id), [2 3]); %#ok<AGROW>
      end
    else
      fin_t = fopen([tform_dir, image_prefixes{tile_id}, ...
        '.global_patchwise_tform.txt'], 'rt');
      transform.transform = fscanf(fin_t, '%g', [6 inf]);
      fclose(fin_t);
      
      n_patch = size(transform.transform, 2);
      for patch_id = 1:n_patch
        if(config.is_verbose)
          fprintf('patch: %d\n', patch_id);
        end
        transform_affine{patch_id} = ...
          reshape(transform.transform(:,patch_id), [3 2])'; %#ok<AGROW>
      end
    end
    
    % add this tile's patches to the section image.
    tile_image = zeros(canvas_height, canvas_width)-1;
    file_name_prefix = [fold_dir, image_prefixes{tile_id}, '.fold_mask'];
    fold_mask = load_fold_mask(file_name_prefix, config);
    for patch_id = 1:size(xdata,2) %#ok<USENS>
      if(isempty(xdata{tile_id, patch_id}))
        continue;
      end
      image = (fold_mask==patch_id)*patch_id;
      image(fold_mask~=patch_id) = -1;
      [image_t, xdata_1, ydata_1] = ...
        imtransform(image, tforms{tile_id, patch_id}, 'FillValue', -1); %#ok<USENS>
      temp_image = zeros(canvas_height, canvas_width)-1;
      temp_image(...
        ydata{tile_id, patch_id}(1):ydata{tile_id, patch_id}(2), ...
        xdata{tile_id, patch_id}(1):xdata{tile_id, patch_id}(2)) = ...
        image_t; %#ok<USENS>
      temp_image(tile_image>=0) = tile_image(tile_image>=0);
      tile_image = temp_image;
      
      transform_affine{patch_id}(1,3) = transform_affine{patch_id}(1,3) - ...
        xdata_1(1) + xdata{tile_id, patch_id}(1); %#ok<AGROW>
      transform_affine{patch_id}(2,3) = transform_affine{patch_id}(2,3) - ...
        ydata_1(1) + ydata{tile_id, patch_id}(1); %#ok<AGROW>
    end
    
    annotations_file_name = [annot_dir, image_prefixes{tile_id}, ...
      config.proofreader.annotations.suffix, '.txt'];
    if(exist(annotations_file_name, 'file')==2)
      fprintf('reading: %s\n', annotations_file_name);
      [annotations_voxel, annotations_body_segment] = ...
        get_tile_voxel_annotations(annotations_file_name);
      
      % Voxel annotations: transform the voxel annotations to global
      % registration coordinates and include them in the section's
      % annotation file if that point in the container patch is visible in
      % the stitched image. 
      if(~isempty(annotations_voxel))
        for a_id = 1:length(annotations_voxel)
          patch_id = fold_mask(annotations_voxel(a_id).y, annotations_voxel(a_id).x);
          if(patch_id<=0)
            continue;
          end
          xt = round(transform_affine{patch_id}(1,:)*...
            [annotations_voxel(a_id).x, annotations_voxel(a_id).y, 1]');
          yt = round(transform_affine{patch_id}(2,:)*...
            [annotations_voxel(a_id).x, annotations_voxel(a_id).y, 1]');
          annotations_voxel(a_id).x = xt;
          annotations_voxel(a_id).y = yt;
          annotations_voxel(a_id).annotation = transform_tbar_postsynaptic_coordinates(...
            annotations_voxel(a_id).annotation, transform_affine{patch_id}, [], [], 0);
          if(section_image(annotations_voxel(a_id).y, annotations_voxel(a_id).x)<=0)
            % a previous tile did not overlap in this area so include this
            % annotation
            annotations_voxel(a_id).x = annotations_voxel(a_id).x - align_roi.xmin + 1;
            annotations_voxel(a_id).y = annotations_voxel(a_id).y - align_roi.ymin + 1;
            annotations_voxel(a_id).annotation = transform_tbar_postsynaptic_coordinates(...
              annotations_voxel(a_id).annotation, ...
              [1, 0, -align_roi.xmin+1; 0, 1, -align_roi.ymin+1], [], [], 0);

            if(isfield(config.proofreader.roi, 'ymin'))
              annotations_voxel(a_id).x = annotations_voxel(a_id).x - ...
                config.proofreader.roi.xmin + 1;
              annotations_voxel(a_id).y = annotations_voxel(a_id).y - ...
                config.proofreader.roi.ymin + 1;
              annotations_voxel(a_id).annotation = transform_tbar_postsynaptic_coordinates(...
                annotations_voxel(a_id).annotation, ...
                [1, 0, -config.proofreader.roi.xmin+1; 0, 1, -config.proofreader.roi.ymin+1], [], [], 0);
            end
            if(isempty(annotations_voxel_section))
              annotations_voxel_section = annotations_voxel(a_id);
            else
              annotations_voxel_section(end+1) = annotations_voxel(a_id); %#ok<AGROW>
            end
          end
        end
      end
      
      % Body segment annotations: apply label_mapping on segment ids.
      if(~isempty(annotations_body_segment))
        for a_id = 1:length(annotations_body_segment)
          annotations_body_segment(a_id).segments = ...
            label_mapping.label_mappings{tile_id}(1+annotations_body_segment(a_id).segments);
          
          if(isempty(annotations_body_segment_section))
            annotations_body_segment_section = annotations_body_segment(a_id);
          else
            annotations_body_segment_section(end+1) = annotations_body_segment(a_id); %#ok<AGROW>
          end
        end
      end
    end
    
    tile_image(section_image>=0) = section_image(section_image>=0);
    section_image = tile_image;
  end
  
  if(~isempty(annotations_voxel_section) || ...
      ~isempty(annotations_body_segment_section))
    fprintf('dumping annotations\n');
    check_for_dir([get_region_dir(config), config.annotations.dir]);
    annotations_file_name = [get_region_dir(config), config.annotations.dir, ...
      'annotations.', num2str(case_id), '.txt'];
    output_section_annotations(annotations_voxel_section, ...
      annotations_body_segment_section, annotations_file_name);
  end
end

if(config.is_verbose)
  fprintf('STOP: generate_annotations_multi_tile_patches\n');
end
return
end
