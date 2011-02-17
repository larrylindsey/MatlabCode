function align_stack_deformable_mesh_tile_pair_in_plane(config)
% align_stack_deformable_mesh_tile_pair_in_plane(config)
% Compute piecewise affine transformations using deformable mesh between
% tiles within a section. May be used for global alignment in the presence
% of folds.
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
% v5  12312008  modified for within section alignment
%

stack_config = config.stack;
dmesh_config = config.align.in_section_tile_pair_align.deformable_mesh;

image_dir = get_stack_dir(config);
dmesh_dir = get_deformable_mesh_dir(config);

if(~isfield(dmesh_config, 'is_verbose'))
  dmesh_config.is_verbose = true;
end
if(~isfield(dmesh_config, 'is_verbose_figures'))
  dmesh_config.is_verbose_figures = false;
end
if(dmesh_config.is_verbose)
  fprintf('START: align_stack_deformable_mesh_tile_pair_in_plane\n');
end

% Generate alignment transforms between pairs of tiles from consequtive
% plane/section
fold_dir = get_fold_dir(config);
norm_x_corr_dir = get_norm_x_corr_dir(config);
for plane = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(plane);
  if(dmesh_config.is_verbose)
    fprintf('case_id: %d\n', case_id);
  end
  
  [image_prefixes, image_sub_dirs, is_to_be_processed] = ...
    get_image_prefixes_subdirs(config, case_id);
  tforms = get_image_tforms(config, case_id);
  sizes = get_image_sizes(config, case_id);
  
  if(dmesh_config.is_verbose)
    fprintf('loading images and fold masks\n');
  end
  images = {};
  fold_masks = {};
  for j = 1:length(image_prefixes)
    file_name=[image_dir, image_prefixes{j}, stack_config.image_suffix];
    images{j} = im2double(imread(file_name)); %#ok<AGROW>
    if(size(images{j},3)>1)
      images{j} = rgb2gray(images{j}); %#ok<AGROW>
    end
    if(max(images{j}(:))>1)
      images{j} = images{j}/max(images{j}(:)); %#ok<AGROW>
    end
    if(isfield(dmesh_config, 'filter_version') && ~isempty(dmesh_config.filter_version))
      images{j} = filter_image(images{j}, dmesh_config.filter_version); %#ok<AGROW>
    end
    images{j} = uint8(255*images{j}); %#ok<AGROW>
    
    file_name_prefix =[fold_dir, image_prefixes{j}, '.fold_mask'];
    fold_masks{j} = load_fold_mask(file_name_prefix, config); %#ok<AGROW>
  end
  % get the default parameters
  params = dmesh_config.params;
  %now get any special ones for this layer
  if(isfield(dmesh_config, 'exceptions') && ...
      ~isempty(dmesh_config.exceptions))
    exc = dmesh_config.exceptions;
    for e=1:size(exc,1)
      if (exc{e,1}(1) == case_id && exc{e,1}(2) == case_id)
        params = [params, ' ', exc{e,2}]; %#ok<AGROW>
      end
    end
  end
  
  %
  % align using deformable mesh if adjacent
  %
  for tile_1 = 1:length(images)
    if(isempty(images{tile_1}))
      continue;
    end
    for tile_2 = tile_1+1:length(images)
      if(isempty(images{tile_2}))
        continue;
      end
    
      if(dmesh_config.is_verbose)
        fprintf('%d,%d ', tile_1, tile_2);
        fprintf('--- Tile pair ---\n');
        fprintf('%s\n%s\n', image_prefixes{tile_1}, ...
          image_prefixes{tile_2});
      end
      
      if(is_to_be_processed(tile_1)==0 && is_to_be_processed(tile_2)==0)
        if(dmesh_config.is_verbose)
          fprintf('both tiles have is_to_be_processed=false, skipping this pair\n');
        end
        continue;
      end

      file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
        image_prefixes{tile_1}, image_prefixes{tile_2}, 'dt.');
      file_name = get_storage_file_name([file_name_suffix, '.mat']);
      if(dmesh_config.is_verbose)
        fprintf('Deleting pre-existing transform file.\n');
      end
      delete(file_name);
      
      if(are_overlapping_affine_tforms(tforms{tile_1}, tforms{tile_2}, ...
          sizes{tile_1}, sizes{tile_2})==0)
        if(dmesh_config.is_verbose)
          fprintf('Tile pair not overlapping, skipping this pair\n');
        end
        continue;
      end
      
      transforms_tp = [];
      
      if(isfield(dmesh_config, 'use_precomputed_overlap') && ...
          dmesh_config.use_precomputed_overlap)
        file_name_suffix = get_file_name_from_tuple(norm_x_corr_dir, ...
          image_prefixes{tile_1}, image_prefixes{tile_2}, 'nc.');
        file_name = [file_name_suffix, '.mat'];
        transform_norm_x_corr = load2(file_name, 'transforms_tp');
        if(isempty(transform_norm_x_corr.transforms_tp))
          continue;
        end
        [overlap_mask_1, overlap_mask_2] = get_overlap_mask(...
          size(images{tile_1}), size(images{tile_2}), ...
          transform_norm_x_corr.transforms_tp);
        overlap_mask_1 = imdilate(overlap_mask_1>0, strel('disk', 10));
        overlap_mask_2 = imdilate(overlap_mask_2>0, strel('disk', 10));
        fold_mask_1 = fold_masks{tile_1} .* uint8(overlap_mask_1);
        fold_mask_2 = fold_masks{tile_2} .* uint8(overlap_mask_2);
      else
        fold_mask_1 = fold_masks{tile_1};
        fold_mask_2 = fold_masks{tile_2};
      end

      [map_mask, transforms] = tile_align_deformable_mesh(...
        images{tile_1}, fold_mask_1, images{tile_2}, fold_mask_2, 1, params);
      
      if(~isempty(transforms))
        % the c-code is row-major where as MATLAB is column-major
        transforms_tp.map_mask = map_mask;
        transforms_tp.transforms = transforms([4 3 2 1 6 5]', :);
        if(dmesh_config.is_verbose_figures)
          display_image = repmat(images{tile_1}, [1 1 3]);
          m1 = [diff(double(map_mask)); zeros(1,size(map_mask,2))];
          m2 = [diff(double(map_mask), [], 2), zeros(size(map_mask,1),1)];
          m = max(abs(m1), abs(m2));
          display_image(:,:,3) = 255*imdilate(m, strel('disk', 10));
          figure;
          subplot(1,2,1); imshow(double(display_image)/255);
          subplot(1,2,2); imshow(images{tile_2});
          title(sprintf('plane %d, tile1 %d, tile2 %d\n', case_ids_ref(plane), tile_1, ...
            tile_2));
        end
      else
        transforms_tp.transforms = [];
        transforms_tp.map_mask = [];
      end
      
      file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
        image_prefixes{tile_1}, image_prefixes{tile_2}, 'dt.');
      file_name = [file_name_suffix, '.mat'];
      save2(file_name, 'transforms_tp');
    end
  end
  
  if(dmesh_config.is_verbose)
    fprintf('\n');
  end
end

if(dmesh_config.is_verbose)
  fprintf('STOP: align_stack_deformable_mesh_tile_pair_in_plane\n');
end

return
end
