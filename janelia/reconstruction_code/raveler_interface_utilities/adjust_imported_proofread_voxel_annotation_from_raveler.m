function adjust_imported_proofread_voxel_annotation_from_raveler(config)
% adjust_imported_proofread_voxel_annotation_from_raveler(config)
% Adjusts imported proofread annotations to account for overlapping region in
% tiles.
%
% Shiv N. Vitaladevuni
% v0  03192010  init. code
%

fprintf('START: adjust_imported_proofread_voxel_annotation_from_raveler\n');
if(strcmp(config.proofreader.method, 'Raveler')~=1)
  error('Proofreader method is not Raveler');
end

stack_config = config.stack;
import_config = config.proofreader.import;
raveler_config = import_config.Raveler;

if(~isfield(import_config, 'is_verbose'))
  import_config.is_verbose = true;
end
if(~isfield(import_config, 'is_verbose_figures'))
  import_config.is_verbose_figures = false;
end

if(import_config.is_verbose)
  fprintf('Adjusting proofread annotations from proofreader in areas of overlapping tiles..\n');
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
  config.annotations.dir, import_config.dir];
check_for_dir(save_dir);
import_suffix = ['.prfd_annot', '.', import_config.proofread_name];
save_suffix = ['.prfrd_adj_annot', '.', import_config.proofread_name];

[junk, annotations_body] = get_proofread_voxel_annotations([...
  raveler_config.proofread_data_dir, raveler_config.annotation_file_name]); %#ok<ASGLU>

fin_seg_2_body = fopen([raveler_config.proofread_data_dir, ...
  raveler_config.segment_to_body_map_file_name], 'rt');
seg_2_body = fscanf(fin_seg_2_body, '%d', [2 inf])';
seg_2_body_map = [];
seg_2_body_map(1+seg_2_body(:,1)) = seg_2_body(:,2);
fclose(fin_seg_2_body);

for layer_id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(layer_id);
  fprintf('plane: %d\n', case_id);
  [image_prefixes, image_sub_dirs] = ...
    get_image_prefixes_subdirs(config, case_id);
  
  fprintf(['Copying imported superpixel maps onto within section', ...
    ' joint alignment ...\n']);
  minx = inf;
  miny = inf;
  maxx = -inf;
  maxy = -inf;
  tforms = {};
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
  
  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))
      continue;
    end
    if(import_config.is_verbose)
      fprintf('tile:%d\n', tile_id);
    end
    
    pid = find([config.region.region_structure.planes(:).z]==case_id);
    sizes{tile_id} = ...
      [config.region.region_structure.planes(pid).tilt_planes(1).tiles(tile_id).height, ...
      config.region.region_structure.planes(pid).tilt_planes(1).tiles(tile_id).width]; %#ok<AGROW>
    
    load2([tform_dir, image_prefixes{tile_id}, image_set_string, ...
      tform_suffix, '.mat'], 'transform_p');
    if(size(transform_p,1)~=6)
      transform_p = transform_p';
    end
    transforms_affine = [transforms_affine, transform_p]; %#ok<AGROW>
    
    tforms{tile_id} = ...
      maketform('affine', reshape(transform_p, [2 3])'); %#ok<AGROW>
    [tiles_t{tile_id}, xdata{tile_id}, ydata{tile_id}] = ...
      imtransform(ones(sizes{tile_id}), tforms{tile_id}, 'nearest', ...
      'FillValues', -1); %#ok<AGROW>
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
  for tile_id = 1:length(tiles_t)
    if(isempty(tiles_t{tile_id}))
      continue;
    end
    xdata{tile_id} = round(xdata{tile_id}-minx+1); %#ok<AGROW>
    ydata{tile_id} = round(ydata{tile_id}-miny+1); %#ok<AGROW>
    
    transforms_affine(5, tile_id) = transforms_affine(5,tile_id) - ...
      xdata_1{tile_id}(1) + xdata{tile_id}(1); %#ok<AGROW>
    transforms_affine(6,tile_id) = transforms_affine(6,tile_id) - ...
      ydata_1{tile_id}(1) + ydata{tile_id}(1); %#ok<AGROW>
  end
  maxx = round(maxx - minx + 1);
  maxy = round(maxy - miny + 1);
  
  % map each tile's annotation into the within-section alignment canvas.
  if(import_config.is_verbose)
    fprintf('Mapping initial imported proofread annotations to within-section alignment canvas... ');
  end
  annotations_voxel = {};
  for tile_id = 1:length(tiles_t)
    if(isempty(tiles_t{tile_id}))
      continue;
    end
    annot_file_name = [save_dir, image_prefixes{tile_id}, import_suffix, '.txt'];
    if(exist(annot_file_name, 'file')==2)
      annotations_voxel{tile_id} = read_imported_proofread_voxel_annotations(...
        annot_file_name); %#ok<AGROW>
      for i = 1:length(annotations_voxel{tile_id})
        T = reshape(transforms_affine(:, tile_id), [2 3]);
        annotations_voxel{tile_id}(i).xt = T(1,:) * ...
          [annotations_voxel{tile_id}(i).x, annotations_voxel{tile_id}(i).y, 1]'; %#ok<AGROW>
        annotations_voxel{tile_id}(i).yt = T(2,:) * ...
          [annotations_voxel{tile_id}(i).x, annotations_voxel{tile_id}(i).y, 1]'; %#ok<AGROW>
        
        % If a T-bar then get the structure information of the T-bar to
        % transform the post-synaptic partners' coordinates.
        annotations_voxel{tile_id}(i).annotation = transform_tbar_postsynaptic_coordinates(...
          annotations_voxel{tile_id}(i).annotation, T, [], [], 0); %#ok<AGROW>
      end
    else
      annotations_voxel{tile_id} = []; %#ok<AGROW>
    end
  end
  
  % For each tile, see if any annotation falls on it. If so, copy it into
  % the tile's adjusted imported annotations.
  if(import_config.is_verbose)
    fprintf('Copy annotations in case of overlap ...\n');
  end
  for tile_id = 1:length(tiles_t)
    if(isempty(tiles_t{tile_id}))
      continue;
    end
    
    % image of the tile in the within-section alignment canvas
    tile_image = zeros(maxy, maxx)-1;
    tile_image(ydata{tile_id}(1):ydata{tile_id}(2), ...
      xdata{tile_id}(1):xdata{tile_id}(2)) = tiles_t{tile_id};
    
    % See if any annotation falls on it. If so, copy it into
    % the tile's adjusted imported annotations.
    annotations_voxel_adj = [];
    for tile_id_1 = 1:length(tiles_t)
      if(isempty(tiles_t{tile_id_1}))
        continue;
      end
      for i = 1:length(annotations_voxel{tile_id_1})
        if(tile_image(round(annotations_voxel{tile_id_1}(i).yt), round(annotations_voxel{tile_id_1}(i).xt))==1)
          if(isempty(annotations_voxel_adj))
            annotations_voxel_adj = annotations_voxel{tile_id_1}(i);
          else
            annotations_voxel_adj(end+1) = annotations_voxel{tile_id_1}(i); %#ok<AGROW>
          end
          T = inv([reshape(transforms_affine(:, tile_id), [2 3]); 0 0 1]);
          
          annotations_voxel_adj(end).x = round(T(1,:) * ...
            [annotations_voxel{tile_id_1}(i).xt, annotations_voxel{tile_id_1}(i).yt, 1]');
          annotations_voxel_adj(end).y = round(T(2,:) * ...
            [annotations_voxel{tile_id_1}(i).xt, annotations_voxel{tile_id_1}(i).yt, 1]');

          % If a T-bar then get the structure information of the T-bar to
          % transform the post-synaptic partners' coordinates.
          annotations_voxel_adj(end).annotation = transform_tbar_postsynaptic_coordinates(...
            annotations_voxel{tile_id_1}(i).annotation, T, [], [], 0);
        end
      end
    end
    
    % Look for any body annotations within this tile
    seg_dir = [get_reconstruction_dir(config), ...
      config.superpixel(1).dir, import_config.dir];
    seg_suffix = ['.prfrd_seg', '.', import_config.proofread_name];
    file_name = [seg_dir, image_prefixes{tile_id}, ...
      seg_suffix, '.mat'];
    seg = load2(file_name, 'superpixel_to_seg_label');
    seg_ids = unique(nonzeros(seg.superpixel_to_seg_label(:,2)));
    segment_body_pairs = [seg_ids, seg_2_body_map(1+seg_ids)'];
    for i = 1:size(segment_body_pairs, 1)
      a_id = find([annotations_body(:).body_id]==segment_body_pairs(i,2),1);
      if(~isempty(a_id))
        segment_body_pairs(i,3) = a_id;
      else
        segment_body_pairs(i,3) = 0;
      end
    end
    
    % Dump the adjusted annotations
    if(~isempty(annotations_voxel_adj) || ...
        (~isempty(segment_body_pairs) && max(segment_body_pairs(:,3))>0))
      fprintf(['Saving: ', save_dir, image_prefixes{tile_id}, ...
        save_suffix, '.txt\n']);
      check_for_dir([save_dir, image_sub_dirs{tile_id}]);
      fout = fopen([save_dir, image_prefixes{tile_id}, ...
        save_suffix, '.txt'], 'wt');
      
      if(~isempty(annotations_voxel_adj))
        for i = 1:length(annotations_voxel_adj)
          fprintf(fout, 'VOXEL %d %d "%s"\n', ...
            annotations_voxel_adj(i).y, ...
            annotations_voxel_adj(i).x, ...
            annotations_voxel_adj(i).annotation);
        end
      end
      
      if(~isempty(segment_body_pairs) && max(segment_body_pairs(:,3))>0)
        body_ids_annotated = unique(...
          segment_body_pairs(segment_body_pairs(:,3)>0,2));
        for b_id = body_ids_annotated'
          p_id = find(segment_body_pairs(:,2)==b_id);
          s_id = unique(segment_body_pairs(p_id,1));
          for i = 1:length(annotations_body(segment_body_pairs(p_id(1),3)).annotations)
            fprintf(fout, 'BODY_SEGMENT {');
            fprintf(fout, '%d ', s_id);
            fprintf(fout, '} "%s"\n', ...
              annotations_body(segment_body_pairs(p_id(1),3)).annotations{i});
          end
        end
      end
      
      fclose(fout);
    end
  end
end

fprintf('STOP: adjust_imported_proofread_voxel_annotation_from_raveler\n');
return
end
