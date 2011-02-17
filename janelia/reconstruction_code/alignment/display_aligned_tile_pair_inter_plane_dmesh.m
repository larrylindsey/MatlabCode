function display_aligned_tile_pair_inter_plane_dmesh(config)
% display_aligned_tile_pair_inter_plane_dmesh(config)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

global config_global

stack_config = config.stack;
dmesh_config = config.align.linkage_align.deformable_mesh;
if(~isfield(config.align.linkage_align, 'display'))
  display_config = [];
else
  display_config = config.align.linkage_align.display;
end

if(dmesh_config.is_verbose)
  fprintf('START: display_aligned_tile_pair_inter_plane_dmesh\n');
end

dmesh_dir = get_deformable_mesh_dir(config);

[images_1, image_prefixes_1] = ...
  get_image_from_stack(config, stack_config.case_ids(1));
for i = 2:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  fprintf('%d:', case_id);
  if(mod(i,20)==0)
    fprintf('\n');
  end;

  [images_2, image_prefixes_2] = ...
    get_image_from_stack(config, case_id);
  
  % display alignment between pairs of tiles
  if(~isfield(display_config, 'tile_1'))
    tiles_1 = 1:length(image_prefixes_1);
  else
    tiles_1 = display_config.tile_1;
  end
  if(~isfield(display_config, 'tile_2'))
    tiles_2 = 1:length(image_prefixes_2);
  else
    tiles_2 = display_config.tile_2;
  end
  for tile_1 = tiles_1
    for tile_2 = tiles_2
      fprintf('--- Tile pair ---\n%d,%d\n', tile_1, tile_2);
      fprintf('%s\n%s\n', image_prefixes_1{tile_1}, image_prefixes_2{tile_2});
      
      fprintf('Loading transform from tile 2 to tile 1\n');
      if(isfield(dmesh_config, 'input_tform_format') && ...
          strcmp(dmesh_config.input_tform_format, 'tif_txt'))
        % Read in piecewise transformations as a TIFF file of the map mask
        % and a text file of the affine transformations.
        fprintf('Trying to load transform txt file\n');
        file_name_prefix = get_file_name_from_tuple(dmesh_dir, ...
          image_prefixes_2{tile_2}, image_prefixes_1{tile_1}, 'dt.');
        file_name_tforms = get_storage_file_name([file_name_prefix, '.tforms.txt']);
        transforms_tp = [];
        if(exist(file_name_tforms, 'file')~=2)
          fprintf('Transform txt file does not exist.\n');
        else
          fprintf('Reading transform txt file.\n');
          fin_tforms = fopen(file_name_tforms, 'rt');
          transforms_tp.transforms = fscanf(fin_tforms, '%g', [6, inf]);
          fclose(fin_tforms);
          if(isempty(transforms_tp.transforms))
            transforms_tp.map_mask = [];
            fprintf('Transform txt file was empty.\n');
          else
            fprintf('Successfully read transform txt file\n');
          end
        end
        fprintf('Trying to load transform tif file\n');
        file_name_map = get_storage_file_name([file_name_prefix, '.map.tif']);
        if(exist(file_name_map, 'file')~=2)
          fprintf('Transform tif file does not exist\n');
        else
          transforms_tp.map_mask = uint16(imread(file_name_map));
          fprintf('Successfully read transform tif file\n');
        end
      else
        % Read in piecewise transformations as a MATLAB .mat file
        file_name_prefix = get_file_name_from_tuple(dmesh_dir, ...
          image_prefixes_2{tile_2}, image_prefixes_1{tile_1}, 'dt.');
        file_name = [file_name_prefix, '.mat'];
        fprintf('Trying to load transform mat file\n');
        transforms_tp = [];
        try
          load2(file_name,'transforms_tp');
          fprintf('Successfully read transform mat file\n');
        catch %#ok<CTCH>
          fprintf('Could not load transform mat file.\n');
        end
      end
      transforms_tp_rev = transforms_tp;

      fprintf('Loading transform from tile 1 to tile 2\n');
      if(isfield(dmesh_config, 'input_tform_format') && ...
          strcmp(dmesh_config.input_tform_format, 'tif_txt'))
        % Read in piecewise transformations as a TIFF file of the map mask
        % and a text file of the affine transformations.
        fprintf('Trying to load transform txt file\n');
        file_name_prefix = get_file_name_from_tuple(dmesh_dir, ...
          image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'dt.');
        file_name_tforms = get_storage_file_name([file_name_prefix, '.tforms.txt']);
        if(exist(file_name_tforms, 'file')~=2)
          fprintf('Transform txt file does not exist.\n');
        else
          fprintf('Reading transform txt file.\n');
          transforms_tp = [];
          fin_tforms = fopen(file_name_tforms, 'rt');
          transforms_tp.transforms = fscanf(fin_tforms, '%g', [6, inf]);
          fclose(fin_tforms);
          if(isempty(transforms_tp.transforms))
            transforms_tp.map_mask = [];
            fprintf('Transform txt file was empty.\n');
          else
            fprintf('Successfully read transform txt file\n');
          end
        end
        fprintf('Trying to load transform tif file\n');
        file_name_map = get_storage_file_name([file_name_prefix, '.map.tif']);
        if(exist(file_name_map, 'file')~=2)
          fprintf('Transform tif file does not exist\n');
        else
          transforms_tp.map_mask = uint16(imread(file_name_map));
          fprintf('Successfully read transform tif file\n');
        end
      else
        % Read in piecewise transformations as a MATLAB .mat file
        file_name_prefix = get_file_name_from_tuple(dmesh_dir, ...
          image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'dt.');
        file_name = [file_name_prefix, '.mat'];
        fprintf('Trying to load transform mat file\n');
        try
          load2(file_name,'transforms_tp');
          fprintf('Successfully read transform mat file\n');
        catch %#ok<CTCH>
          fprintf('Could not load transform mat file. Going to next pair\n');
        end
      end
      
      [coords_x_1, coords_y_1] = meshgrid(1:size(images_1{tile_1},2), ...
        1:size(images_1{tile_1},1));
      [coords_x_2, coords_y_2] = meshgrid(1:size(images_1{tile_2},2), ...
        1:size(images_1{tile_2},1));

      if(~isempty(transforms_tp))
        transforms_tp.map_mask(transforms_tp.map_mask>0) = ...
          transforms_tp.map_mask(transforms_tp.map_mask>0) - ...
          config_global.TRANSFORMATION_ID_OFFSET + 1;
        t1 = [0 transforms_tp.transforms(1,:)];
        t2 = [0 transforms_tp.transforms(2,:)];
        t3 = [0 transforms_tp.transforms(3,:)];
        t4 = [0 transforms_tp.transforms(4,:)];
        t5 = [0 transforms_tp.transforms(5,:)];
        t6 = [0 transforms_tp.transforms(6,:)];
        coords_x_1_t = t1(transforms_tp.map_mask+1).*coords_x_1 + ...
          t3(transforms_tp.map_mask+1).*coords_y_1 + t5(transforms_tp.map_mask+1);
        coords_x_1_t = round(coords_x_1_t);
        coords_x_1_t(coords_x_1_t<1) = 1;
        coords_x_1_t(coords_x_1_t>size(images_2{tile_2},2)) = size(images_2{tile_2},2);
        coords_y_1_t = t2(transforms_tp.map_mask+1).*coords_x_1 + ...
          t4(transforms_tp.map_mask+1).*coords_y_1 + t6(transforms_tp.map_mask+1);
        coords_y_1_t = round(coords_y_1_t);
        coords_y_1_t(coords_y_1_t<1) = 1;
        coords_y_1_t(coords_y_1_t>size(images_2{tile_2},1)) = size(images_2{tile_2},1);
        d = repmat(1-images_1{tile_1}, [1 1 3]);
        d1 = images_2{tile_2}(sub2ind(size(images_1{tile_1}), ...
          coords_y_1_t, coords_x_1_t));
        d(:,:,2) = 1-d1;
        figure;
        imshow(imresize(d, 1/2));
        title(['tile1 to tile2 ', num2str(tile_1), ' ', image_prefixes_1{tile_1}(1:30), ...
          num2str(tile_2), ' ', image_prefixes_2{tile_2}(1:30)]);
      end
      
      if(~isempty(transforms_tp_rev))
        transforms_tp_rev.map_mask(transforms_tp_rev.map_mask>0) = ...
          transforms_tp_rev.map_mask(transforms_tp_rev.map_mask>0) - ...
          config_global.TRANSFORMATION_ID_OFFSET + 1;
        t1 = [0 transforms_tp_rev.transforms(1,:)];
        t2 = [0 transforms_tp_rev.transforms(2,:)];
        t3 = [0 transforms_tp_rev.transforms(3,:)];
        t4 = [0 transforms_tp_rev.transforms(4,:)];
        t5 = [0 transforms_tp_rev.transforms(5,:)];
        t6 = [0 transforms_tp_rev.transforms(6,:)];
        coords_x_2_t = t1(transforms_tp_rev.map_mask+1).*coords_x_2 + ...
          t3(transforms_tp_rev.map_mask+1).*coords_y_2 + t5(transforms_tp_rev.map_mask+1);
        coords_x_2_t = round(coords_x_2_t);
        coords_x_2_t(coords_x_2_t<1) = 1;
        coords_x_2_t(coords_x_2_t>size(images_1{tile_1},2)) = size(images_1{tile_1},2);
        coords_y_2_t = t2(transforms_tp_rev.map_mask+1).*coords_x_2 + ...
          t4(transforms_tp_rev.map_mask+1).*coords_y_2 + t6(transforms_tp_rev.map_mask+1);
        coords_y_2_t = round(coords_y_2_t);
        coords_y_2_t(coords_y_2_t<1) = 1;
        coords_y_2_t(coords_y_2_t>size(images_1{tile_1},1)) = size(images_1{tile_1},1);
        d = 1-repmat(images_2{tile_2}, [1 1 3]);
        d2 = images_1{tile_1}(sub2ind(size(images_2{tile_2}), ...
          coords_y_2_t, coords_x_2_t));
        d(:,:,2) = 1-d2;
        figure;
        imshow(imresize(d, 1/2));
        title(['tile2 to tile1 ', num2str(tile_2), ' ', image_prefixes_2{tile_2}(1:30), ...
          num2str(tile_1), ' ', image_prefixes_1{tile_1}(1:30)]);
      end
    end
  end

  images_1 = images_2;  
  image_prefixes_1 = image_prefixes_2;
  fprintf('\n');
end;

if(dmesh_config.is_verbose)
  fprintf('STOP: display_aligned_tile_pair_inter_plane_dmesh\n');
end
