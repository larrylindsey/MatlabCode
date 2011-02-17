function adjust_imported_proofread_superpixel_map_tile_patches(config)
% adjust_imported_proofread_superpixel_map(config)
% Adjusts imported suerpixel maps to account for overlapping region in
% tiles.
% Label '0' is written on all undefined pixels.
%
% Shiv N. Vitaladevuni
% v0  04062009  init. code
%

fprintf('START: adjust_imported_proofread_superpixel_map_tile_patches\n');
if(strcmp(config.proofreader.method, 'Raveler')~=1)
  error('Proofreader method is not Raveler');
end

fprintf('Importing proofread data from proofreader ..\n');

stack_config = config.stack;
import_config = config.proofreader.import;
if(~isfield(import_config, 'is_verbose'))
  import_config.is_verbose = true;
end
if(~isfield(import_config, 'is_verbose_figures'))
  import_config.is_verbose_figures = false;
end

switch(config.align.in_section_align.method)
  case 'SIFT'
    tform_dir = get_sift_dir(config);
    tform_suffix = '.sift_transform.in_plane';
  case 'deformable_mesh'
    tform_dir = get_deformable_mesh_dir(config);
    tform_suffix = '.dmesh_affine.in_plane';
end

save_dir = [get_reconstruction_dir(config), ...
  config.superpixel(1).dir, import_config.dir];
check_for_dir(save_dir);
import_suffix = ['.prfrd_csp', '.', import_config.proofread_name];
save_suffix = ['.prfrd_asp', '.', import_config.proofread_name];
for layer_id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(layer_id);
  fprintf('plane: %d\n', case_id);
  [images, image_prefixes, image_sub_dirs] = ...
    get_image_from_stack(config, case_id);

  fprintf(['Copying imported superpixel maps onto within section', ...
    ' joint alignment ...\n']);
  minx = inf;
  miny = inf;
  maxx = -inf;
  maxy = -inf;
  tforms = {};
  segs_t = {};
  tiles_t = {};
  xdata = {};
  ydata = {};
  sizes = {};
  transforms_affine = [];
  if(import_config.is_verbose)
    fprintf('Computing transforms ...\n');
  end
  
  % Generate a substring for the file name that accounts for the images
  % involved in the joint alignment
  image_set_string = get_set_string(image_prefixes);
  
  for tile_id = 1:length(images)
    if(isempty(images{tile_id}))
      continue;
    end
    if(import_config.is_verbose)
      fprintf('tile:%d\n', tile_id);
    end
    file_name = [save_dir, image_prefixes{tile_id}, ...
      import_suffix, '.mat'];
    seg = load2(file_name, 'label_map');
    sizes{tile_id} = size(seg.label_map); %#ok<AGROW>
    load2([tform_dir, image_prefixes{tile_id}, image_set_string, ...
      tform_suffix, '.mat'], 'transform_p');
    if(size(transform_p,1)~=6)
      transform_p = transform_p';
    end
    transforms_affine = [transforms_affine, transform_p]; %#ok<AGROW>

    tforms{tile_id} = ...
      maketform('affine', reshape(transform_p, [2 3])'); %#ok<AGROW>
    temp = tile_id*(seg.label_map>=0);
    temp(temp==0) = -1;
    temp = imtransform(temp, tforms{tile_id}, 'nearest', ...
      'FillValues', -1);
    [segs_t{tile_id}, xdata{tile_id}, ydata{tile_id}] = ...
      imtransform(double(seg.label_map), tforms{tile_id}, 'nearest', ...
      'FillValues', -1); %#ok<AGROW>
    segs_t{tile_id}(temp<=0) = -1; %#ok<AGROW>
    tiles_t{tile_id} = temp; %#ok<AGROW>
    minx = min([minx, xdata{tile_id}]);
    miny = min([miny, ydata{tile_id}]);
    maxx = max([maxx, xdata{tile_id}]);
    maxy = max([maxy, ydata{tile_id}]);
  end
  if(import_config.is_verbose)
    fprintf('done\n');
  end

  xdata_1 = xdata;
  ydata_1 = ydata;
  for tile_id = 1:length(segs_t)
    if(isempty(segs_t{tile_id}))
      continue;
    end
    xdata{tile_id} = round(xdata{tile_id}-minx+1); %#ok<AGROW>
    ydata{tile_id} = round(ydata{tile_id}-miny+1); %#ok<AGROW>
  end
  maxx = round(maxx - minx + 1);
  maxy = round(maxy - miny + 1);

  if(import_config.is_verbose)
    fprintf('Rendering ... ');
  end
  stiched_sp_inplane = zeros(maxy, maxx)-1;
  tile_image = zeros(maxy, maxx)-1;
  for tile_id = 1:length(segs_t)
    if(isempty(segs_t{tile_id}))
      continue;
    end
    temp_image = zeros(maxy, maxx)-1;
    temp_image(ydata{tile_id}(1):ydata{tile_id}(2), ...
      xdata{tile_id}(1):xdata{tile_id}(2)) = segs_t{tile_id};
    temp_image(stiched_sp_inplane>=0) = ...
      stiched_sp_inplane(stiched_sp_inplane>=0);
    stiched_sp_inplane = temp_image;

    temp_image = zeros(maxy, maxx)-1;
    temp_image(ydata{tile_id}(1):ydata{tile_id}(2), ...
      xdata{tile_id}(1):xdata{tile_id}(2)) = tiles_t{tile_id};
    temp_image(tile_image>=0) = ...
      tile_image(tile_image>=0);
    tile_image = temp_image;
    
    transforms_affine(5, tile_id) = transforms_affine(5,tile_id) - ...
      xdata_1{tile_id}(1) + xdata{tile_id}(1); %#ok<AGROW>
    transforms_affine(6,tile_id) = transforms_affine(6,tile_id) - ...
      ydata_1{tile_id}(1) + ydata{tile_id}(1); %#ok<AGROW>
  end
  stiched_sp_inplane(stiched_sp_inplane<0) = 0;
  tile_image(tile_image<0) = 0;
  if(import_config.is_verbose_figures)
    figure(layer_id*10);
    imshow(stiched_sp_inplane, []);
    title(['Rendered aligned imported superpixel maps. Layer: ', ...
      num2str(case_id)]);
    figure(layer_id*10+1);
    imshow(tile_image, []);
    title(['Rendered aligned tiles. Layer: ', num2str(case_id)]);
  end
  if(import_config.is_verbose)
    fprintf('done\n');
  end
  
  for tile_id = 1:length(segs_t)
    if(isempty(segs_t{tile_id}))
      continue;
    end
    if(import_config.is_verbose)
      fprintf('tile: %d\n', tile_id);
    end
    fprintf(['Copying overlapping region''s superpixels to images', ...
      ' underneath.\n']);
    overlap_sp_map = get_label_piecewise_affine_transformation_mask(...
      double(stiched_sp_inplane), uint32(tile_id*ones(sizes{tile_id})), ...
      double(transforms_affine), -1);
    overlap_map = get_label_piecewise_affine_transformation_mask(...
      double(tile_image), uint32(tile_id*ones(sizes{tile_id})), ...
      double(transforms_affine), -1);
    overlap_map(overlap_map==tile_id) = 0;
    overlap_sp_map = int32(overlap_sp_map);
    overlap_sp_map(overlap_map<=0) = -1;
    if(import_config.is_verbose_figures)
      figure(layer_id*10+1+tile_id);
      imshow(images{tile_id});
      [py, px] = find(overlap_sp_map==0);
      hold on; plot(px, py, '.'); hold off;
      title(['Superpixel map in the overlapping area. Layer: ', ...
        num2str(case_id), ' tile: ', num2str(tile_id)]);
    end
    fprintf('done.\n');

    fprintf('Adjusting superpixels in overlapping regions.\n');
    image = im2double(images{tile_id});
    image = image/max(image(:));
    image = medfilt2(image, [3 3]);
    boundary = 1 - image;
    boundary(overlap_map<=0) = 1;
    seed_mask = overlap_sp_map;
    seed_mask(seed_mask<0) = 0;
    seed_mask = imerode(seed_mask, strel('disk', 3));
    adjusted_overlap_sp_map = watershed_seeded(boundary, seed_mask);

    valid_overlap_mask = bwareaopen(imclose(overlap_sp_map>0, strel('disk', 4)), ...
      5000);
    adjusted_overlap_sp_map = valid_overlap_mask.*adjusted_overlap_sp_map;
    
    file_name = [save_dir, image_prefixes{tile_id}, ...
      import_suffix, '.mat'];
    seg = load2(file_name, 'label_map');
    label_map = adjusted_overlap_sp_map;
    label_map(overlap_sp_map==-1) = ...
      seg.label_map(overlap_sp_map==-1);
    if(import_config.is_verbose_figures)
      figure(layer_id*20+1+tile_id);
      imshow(images{tile_id});
      [py, px] = find(label_map==0);
      hold on; plot(px, py, '.'); hold off;
      title(['Adjusted superpixel map. Layer: ', ...
        num2str(case_id), ' tile: ', num2str(tile_id)]);
    end
    
    label_map(label_map<0) = 0;
    label_mask = label_map>0;
    label_mask = bwareaopen(label_mask, 10);
    label_map(label_mask==0) = 0;
    fprintf('done.\n');
    
    fprintf('Saving.\n');
    check_for_dir([save_dir, image_sub_dirs{tile_id}]);
    file_name = [save_dir, image_prefixes{tile_id}, ...
      save_suffix, '.mat'];
    locked_labels = unique(label_map(label_map>0)); %#ok<NASGU>
    save2(file_name, 'label_map', 'locked_labels');
    fprintf('done.\n');
  end
  
  
  fprintf('done.\n');
end

fprintf('STOP: adjust_imported_proofread_superpixel_map_tile_patches\n');
return
end
