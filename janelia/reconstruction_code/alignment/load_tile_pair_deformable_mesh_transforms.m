function [transforms_tp, transforms_tp_rev] = ...
  load_tile_pair_deformable_mesh_transforms(image_prefix_1, image_prefix_2, ...
  dmesh_config, dmesh_dir)

if(isfield(dmesh_config, 'input_tform_format') && ...
    strcmp(dmesh_config.input_tform_format, 'tif_txt'))
  % Read in piecewise transformations as a TIFF file of the map mask
  % and a text file of the affine transformations.
  fprintf('Trying to load transform txt file\n');
  file_name_prefix = get_file_name_from_tuple(dmesh_dir, ...
    image_prefix_2, image_prefix_1, 'dt.');
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
    image_prefix_2, image_prefix_1, 'dt.');
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
    image_prefix_1, image_prefix_2, 'dt.');
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
    image_prefix_1, image_prefix_2, 'dt.');
  file_name = [file_name_prefix, '.mat'];
  fprintf('Trying to load transform mat file\n');
  try
    load2(file_name,'transforms_tp');
    fprintf('Successfully read transform mat file\n');
  catch %#ok<CTCH>
    fprintf('Could not load transform mat file. Going to next pair\n');
  end
end
return
end
