function align_stack_deformable_mesh_tile_pair_inter_plane(config)
% align_stack_deformable_mesh_tile_pair_inter_plane(config)
% Compute piecewise affine transformations using deformable mesh between
% pairs of tiles of consequtive sections.
%
% Input:
%   config    config datastructure of the reconstruction
%
% Deformable mesh code in code/aligment/tile_align_fold_deformable_mesh/
% Lou Scheffer
% Visiting Scientist, Janelia Farm Research Campus, HHMI.
%
% Wrapped for pipeline by Shiv Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  08252008  init. code
% v1  08312008  modifications for pipeline - Shiv N. Vitaladevuni, JFRC, HHMI
% v2  09192008  split into SIFT feature point, matching and transform.
% v3  09262008  code modified from in plane alignment to inter plane
%                 tile-pair alignment. 
% v4  12112008  code modified for deformable mesh algorithm
%

stack_config = config.stack;
dmesh_config = config.align.linkage_align.deformable_mesh;

image_dir = get_stack_dir(config);
dmesh_dir = get_deformable_mesh_dir(config);

if(~isfield(dmesh_config, 'is_verbose'))
  dmesh_config.is_verbose = true;
end
if(dmesh_config.is_verbose)
  fprintf('START: align_stack_deformable_mesh_tile_pair_inter_plane\n');
end

if(isfield(dmesh_config, 'use_segmentation_map') && ...
    dmesh_config.use_segmentation_map)
  seg_dir = [get_reconstruction_dir(config), config.segmentation_2D.dir, ...
    dmesh_config.segmentation_method, '/'];
else
  dmesh_config.use_segmentation_map = false;
end

% Generate alignment transforms between pairs of tiles from consequtive
% plane/section
fold_dir = get_fold_dir(config);
norm_x_corr_dir = get_norm_x_corr_dir(config);
images_1 = {};
for plane_1 = 1:length(stack_config.case_ids)-1 % for each plane/section
  case_id_1 = stack_config.case_ids(plane_1);
  if(dmesh_config.is_verbose)
    fprintf('case_id_1: %d\n', case_id_1);
  end

  case_id_2 = case_id_1 + 1;
  if(~ismember(case_id_2, stack_config.case_ids))
    images_1 = {};
    continue;
  end
  
  %
  % load images in first plane if needed
  %
  if(isempty(images_1))
    [image_prefixes_1, image_sub_dirs_1, is_to_be_processed_1] = ...
      get_image_prefixes_subdirs(config, case_id_1);
    tforms_1 = get_image_tforms(config, case_id_1);
    sizes_1 = get_image_sizes(config, case_id_1);
    
    fold_masks_1 = {};
    for j = 1:length(image_prefixes_1)
      if(isempty(image_prefixes_1{j}))
        continue;
      end
      if(~dmesh_config.use_segmentation_map)
        file_name=[image_dir, image_prefixes_1{j}, stack_config.image_suffix];
        images_1{j} = im2double(imread(file_name)); %#ok<AGROW>
        if(size(images_1{j},3)>1)
          images_1{j} = rgb2gray(images_1{j}); %#ok<AGROW>
        end
        if(max(images_1{j}(:))>1)
          images_1{j} = images_1{j}/max(images_1{j}(:)); %#ok<AGROW>
        end
        if(isfield(dmesh_config, 'filter_version') && ~isempty(dmesh_config.filter_version))
          images_1{j} = filter_image(images_1{j}, dmesh_config.filter_version); %#ok<AGROW>
        end
        images_1{j} = uint8(255*images_1{j}); %#ok<AGROW>
      else
        file_name=[seg_dir, image_prefixes_1{j}, ...
          dmesh_config.segmentation_suffix, '.mat'];
        load2(file_name);
        images_1{j} = double(label_map==0); %#ok<AGROW>
        if(isfield(dmesh_config, 'filter_version') && ~isempty(dmesh_config.filter_version))
          images_1{j} = filter_image(images_1{j}, dmesh_config.filter_version); %#ok<AGROW>
        end
      end
      file_name_prefix =[fold_dir, image_prefixes_1{j}, '.fold_mask'];
      fold_masks_1{j} = load_fold_mask(file_name_prefix, config); %#ok<AGROW>
    end
  end
  
  [image_prefixes_2, image_sub_dirs_2, is_to_be_processed_2] = ...
    get_image_prefixes_subdirs(config, case_id_2);
  tforms_2 = get_image_tforms(config, case_id_2);
  sizes_2 = get_image_sizes(config, case_id_2);
  images_2 = {};
  fold_masks_2 = {};
  for j = 1:length(image_prefixes_2)
    if(isempty(image_prefixes_2{j}))
      continue;
    end
    if(~dmesh_config.use_segmentation_map)
      file_name=[image_dir, image_prefixes_2{j}, stack_config.image_suffix];
      images_2{j} = im2double(imread(file_name)); %#ok<AGROW>
      if(size(images_2{j},3)>1)
        images_2{j} = rgb2gray(images_2{j}); %#ok<AGROW>
      end
      if(max(images_2{j}(:))>1)
        images_2{j} = images_2{j}/max(images_2{j}(:)); %#ok<AGROW>
      end
      if(isfield(dmesh_config, 'filter_version') && ~isempty(dmesh_config.filter_version))
        images_2{j} = filter_image(images_2{j}, dmesh_config.filter_version); %#ok<AGROW>
      end
      images_2{j} = uint8(255*images_2{j}); %#ok<AGROW>
    else
      file_name=[seg_dir, image_prefixes_2{j}, ...
        dmesh_config.segmentation_suffix, '.mat'];
      load2(file_name);
      images_2{j} = double(label_map==0); %#ok<AGROW>
      if(isfield(dmesh_config, 'filter_version') && ~isempty(dmesh_config.filter_version))
        images_2{j} = filter_image(images_2{j}, dmesh_config.filter_version); %#ok<AGROW>
      end
    end
    file_name_prefix =[fold_dir, image_prefixes_2{j}, '.fold_mask'];
    fold_masks_2{j} = load_fold_mask(file_name_prefix, config); %#ok<AGROW>
  end
  
  % get the default parameters
  % params1 will be used for the first of the two calls, params2 for the
  % second.  First they both get the general parameters
  params1 = dmesh_config.params;
  params2 = params1;
  %now get any special ones for this layer pair
  if(isfield(dmesh_config, 'exceptions') && ...
      ~isempty(dmesh_config.exceptions))
    exc = dmesh_config.exceptions;
    for e=1:size(exc,1)
      if (exc{e,1}(1) == case_id_1 && exc{e,1}(2) == case_id_2)
        params1 = [params1, ' ', exc{e,2}]; %#ok<AGROW>
      end
      if (exc{e,1}(1) == case_id_2 && exc{e,1}(2) == case_id_1)
        params2 = [params2, ' ', exc{e,2}]; %#ok<AGROW>
      end
    end
  end
  
  %
  % align using deformable mesh if adjacent
  %
  for tile_1 = 1:length(images_1)
    if(isempty(images_1{tile_1}))
      continue;
    end
    for tile_2 = 1:length(images_2)
      if(isempty(images_2{tile_2}))
        continue;
      end
      
      if(dmesh_config.is_verbose)
        fprintf('--- Tile pair ---\n');
        fprintf('%d,%d\n', tile_1, tile_2);
        fprintf('%s\n%s\n', image_prefixes_1{tile_1}, ...
          image_prefixes_2{tile_2});
      end
      
      if(is_to_be_processed_1(tile_1)==0 && is_to_be_processed_2(tile_2)==0)
        if(dmesh_config.is_verbose)
          fprintf('both tiles have is_to_be_processed=false, skipping this pair\n');
        end
        continue;
      end
      
      if(dmesh_config.is_verbose)
        fprintf('Deleting pre-existing transform file.\n');
      end
      file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
        image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'dt.');
      file_name = get_storage_file_name([file_name_suffix, '.mat']);
      delete(file_name);
      file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
        image_prefixes_2{tile_2}, image_prefixes_1{tile_1}, 'dt.');
      file_name = get_storage_file_name([file_name_suffix, '.mat']);
      delete(file_name);
      
      if(are_overlapping_affine_tforms(tforms_1{tile_1}, tforms_2{tile_2}, ...
          sizes_1{tile_1}, sizes_2{tile_2})==0)
        if(dmesh_config.is_verbose)
          fprintf('tiles are not overlapping, skipping this pair\n');
        end
        continue;
      end
      
      if(isfield(dmesh_config, 'use_precomputed_overlap') && ...
          dmesh_config.use_precomputed_overlap)
        file_name_suffix = get_file_name_from_tuple(norm_x_corr_dir, ...
          image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'nc.');
        file_name = [file_name_suffix, '.mat'];
        transform_norm_x_corr = load2(file_name,'transforms_tp');
        if(isempty(transform_norm_x_corr.transforms_tp))
          continue;
        end
        [overlap_mask_1, overlap_mask_2] = get_overlap_mask(...
          size(images_1{tile_1}), size(images_2{tile_2}), ...
          transform_norm_x_corr.transforms_tp);
        overlap_mask_1 = imdilate(overlap_mask_1>0, strel('disk', 10));
        overlap_mask_2 = imdilate(overlap_mask_2>0, strel('disk', 10));
        fold_mask_1 = fold_masks_1{tile_1} .* uint8(overlap_mask_1);
        fold_mask_2 = fold_masks_2{tile_2} .* uint8(overlap_mask_2);
      else
        fold_mask_1 = fold_masks_1{tile_1};
        fold_mask_2 = fold_masks_2{tile_2};
      end
      
      % flip center angle for deformable mesh program if the sections are
      % rotated by 180 degrees.  Since this should be the same for all
      % tiles, this could be done up front
      if(tforms_1{tile_1}(1)*tforms_2{tile_2}(1)>=0)
         p1 = params1;
         p2 = params2;
      else
         p1 = [params1,' -CTR=180'];
         p2 = [params2,' -CTR=180'];
      end
      display(p1);
      display(p2);
      
      % Get transformations in both directions
      tic
      [map_mask, transforms] = tile_align_deformable_mesh(...
        images_1{tile_1}, fold_mask_1, images_2{tile_2}, fold_mask_2, 0, p1);
      toc
      transforms_tp = [];
      transforms_tp.map_mask = map_mask;
      if(~isempty(transforms))
        % the c-code is row-major where as MATLAB is column-major
        transforms_tp.transforms = transforms([4 3 2 1 6 5]', :);
      else
        transforms_tp.transforms = [];
      end
      file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
        image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'dt.');
      file_name = [file_name_suffix, '.mat'];
      save2(file_name, 'transforms_tp');
      
      tic
      [map_mask, transforms] = tile_align_deformable_mesh(...
        images_2{tile_2}, fold_mask_2, images_1{tile_1}, fold_mask_1, 0, p2);
      toc
      transforms_tp = [];
      transforms_tp.map_mask = map_mask;
      if(~isempty(transforms))
        % the c-code is row-major where as MATLAB is column-major
        transforms_tp.transforms = transforms([4 3 2 1 6 5]', :);
      else
        transforms_tp.transforms = [];
      end
      file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
        image_prefixes_2{tile_2}, image_prefixes_1{tile_1}, 'dt.');
      file_name = [file_name_suffix, '.mat'];
      save2(file_name, 'transforms_tp');
    end
  end
  
  if(dmesh_config.is_verbose)
    fprintf('\n');
  end
  
  images_1 = images_2;
  fold_masks_1 = fold_masks_2;
  image_prefixes_1 = image_prefixes_2;
  image_sub_dirs_1 = image_sub_dirs_2; %#ok<NASGU>
  is_to_be_processed_1 = is_to_be_processed_2;
  tforms_1 = tforms_2;
  sizes_1 = sizes_2;
end

if(dmesh_config.is_verbose)
  fprintf('STOP: align_stack_deformable_mesh_tile_pair_inter_plane\n');
end
return
end
