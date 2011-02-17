function align_stack_dmesh_affine_tile_patches(config)
%  align_stack_dmesh_affine_tile_patches(config)
% Compute affine transformations based on correspondences from deformable
% mesh. Break the tile into pieces, each having an affine transformation.
% May be useful in case of folds, tears, etc.
%
% Input:
%   config    config datastructure of the reconstruction
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  010709  init. code
%

stack_config = config.stack;
dmesh_config = config.align.global_align.deformable_mesh;

image_dir = get_stack_dir(config);
dmesh_dir = get_deformable_mesh_dir(config);

if( ~isfield(stack_config, 'image_structure') || isempty(stack_config.image_structure) )
  error('xml file from TrakEM not specified. Could not load tile images and initial transformations.');
end;

if(~isfield(dmesh_config, 'is_verbose'))
  dmesh_config.is_verbose = true;
end
if(~isfield(dmesh_config, 'is_verbose_figures'))
  dmesh_config.is_verbose_figures = false;
end
if(dmesh_config.is_verbose)
  fprintf('Global alignment based on deformable mesh correspondences\n');
end

if(~isfield(dmesh_config, 'min_n_matches'))
  dmesh_config.min_n_matches = 4;
end

s = xmlread([image_dir, stack_config.image_structure]);
xmlwrite('~tmp.xml', s);
xmlstr = fileread('~tmp.xml');
s = xml_parseany(xmlstr);
s_t2_layer = s.t2_layer_set{1}.t2_layer;

if(dmesh_config.is_verbose)
  fprintf('Collecting stack images and initial trakEM transformations ..');
end
if(dmesh_config.is_verbose_figures)
  images = {};
end
image_prefixes = {};
image_sub_dirs = {};
tforms={};
sizes={};
case_ids_ref = -ones(1,length(s_t2_layer)); % record of case_id of each plane to be used when saving
for i = 1:length(s_t2_layer)
  case_id=str2double(s_t2_layer{i}.ATTRIBUTE.z);
  if(ismember(case_id,stack_config.case_ids) || isempty(stack_config.case_ids))
    case_ids_ref(i) = case_id;
    for j = 1:length(s_t2_layer{i}.t2_patch)
      % remove the directory and extension from the file_path and store
      % in image_prefix
      dir_length = strfind(s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path, '/');
      dir_length = [zeros(1, max(0, 3-length(dir_length))), dir_length]; %#ok<AGROW>
      
      % File path in .xml file is
      % <stack_directory>/<grid_id>/<section_id>/<image_file_name>
      if(dmesh_config.is_verbose_figures)
        images{i,j} = im2double(imread(...
          [image_dir, s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path(...
          dir_length(end-2)+1:end)]));
        if(max(images{i,j}(:))>1)
          images{i,j} = images{i,j}/max(images{i,j}(:));
        end
      end
      image_prefixes{i,j} = s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path(dir_length(end-2)+1:end-4);
      image_sub_dirs{i,j} = s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path(dir_length(end-2)+1:dir_length(end));
      
      %Extract transforms
      transform_str = s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.transform;
      if(isempty(transform_str))
        tforms{i,j}=[];
      else
        matrix_str_offset = strfind(s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.transform, '(') + 1;
        matrix_str_end = strfind(s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.transform, ')') - 1;
        xform_matrix = str2num(transform_str(matrix_str_offset:matrix_str_end)); %#ok<ST2NM>
        tforms{i,j} = reshape(xform_matrix, [2 3]);
        
        width_str = s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.width;
        height_str = s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.height;
        sizes{i,j}=[str2double(width_str),str2double(height_str)];
      end
    end
  end
end
if(dmesh_config.is_verbose)
  fprintf('done.\n');
end


%Generate adjacency matrix from trakEM2 data
if(dmesh_config.trakEM2.getAdjacency && ~isempty(tforms))
  if(dmesh_config.is_verbose)
    fprintf('Computing tile adjacency matrix ..\n');
  end
  adjacencies=getSIFTalpha(tforms,sizes);
  if(dmesh_config.is_verbose)
    fprintf('done.\n');
  end
else
  adjacencies=[];
end

if(dmesh_config.save)
  if(~isempty(adjacencies))
    file_name=[dmesh_dir, 'adjacencies.mat'];
    save2(file_name, 'adjacencies');
  end
end

%%%%%%%% Load correspondences from deformable mesh alignments %%%%%%
% % Each tile image is broken into a set of patches according to the fold
% % passing through it. The stack is a union of these sets of patches.
% % An affine transformation is computed for patch in the global reference
% % frame.
% Each patch has a within-tile and global id. This mapping is maintained in
% patch_back_ref{plane,tile} =
%     [global_id_of_patch_1, global_id_of_patch_2, ... ]
if(dmesh_config.is_verbose)
  fprintf('Initializing structures .. \n');
end
n_patch_global = 0;
patch_back_ref = {};
n_patch = [];
fold_dir = get_fold_dir(config);
for plane = 1:size(tforms,1) % for each plane/section
  case_id = case_ids_ref(plane);
  if(case_id<0)
    continue;
  end
  if(dmesh_config.is_verbose)
    fprintf('plane: %d\n', case_id);
  end
  % load fold masks
  fold_masks = struct('fold_mask', []);
  for tile = 1:size(tforms,2)
    if(isempty(tforms{plane,tile}))
      continue;
    end
    file_name_prefix = [fold_dir, image_prefixes{plane, tile}, '.fold_mask'];
    fold_masks(tile).fold_mask = load_fold_mask(file_name_prefix, config);
    n_patch(plane, tile) = double(max(fold_masks(tile).fold_mask(:)));
    patch_back_ref{plane, tile} = ...
      n_patch_global+1 : n_patch_global+n_patch(plane, tile);
    n_patch_global = n_patch_global + n_patch(plane, tile);
  end
end
print_debug('n_patch\n');
disp_debug(n_patch);
if(dmesh_config.is_verbose)
  fprintf('done\n');
end

match_points_a{n_patch_global, n_patch_global} = [];
match_points_b{n_patch_global, n_patch_global} = [];
for plane_1 = 1:size(tforms,1) % for each plane/section
  case_id = case_ids_ref(plane_1);
  if(case_id<0)
    continue;
  end
  if(dmesh_config.is_verbose)
    fprintf('plane1: %d\n', case_id);
  end
  % load fold masks
  fold_masks_1 = struct('fold_mask', []);
  for tile = 1:size(tforms,2)
    if(isempty(tforms{plane_1,tile}))
      continue;
    end
    file_name_prefix = [fold_dir, image_prefixes{plane_1, tile}, '.fold_mask'];
    fold_masks_1(tile).fold_mask = load_fold_mask(file_name_prefix, config);
  end
  
  %
  % deformable mesh transforms for within this plane
  %
  % generate correspondences and split according to folds
  print_debug('Correspondences within section\n');
  for tile_a = 1:size(tforms,2)
    if(isempty(tforms{plane_1,tile_a}))
      continue;
    end
    for tile_b = tile_a+1:size(tforms,2)
      if(isempty(tforms{plane_1,tile_b}))
        continue;
      end
      print_debug('tile pair: %d, %d\n', tile_a, tile_b);
      transforms_tp = [];
      if(isfield(dmesh_config, 'input_tform_format') && ...
          strcmp(dmesh_config.input_tform_format, 'tif_txt'))
        % Read in piecewise transformations as a TIFF file of the map mask
        % and a text file of the affine transformations.
        file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
          image_prefixes{plane_1, tile_a}, image_prefixes{plane_1, tile_b}, 'dt.');
        file_name_tforms = get_storage_file_name([file_name_suffix, '.tforms.txt']);
        if(exist(file_name_tforms, 'file')==2)
          transforms_tp = [];
          fprintf('%s\n', file_name_tforms);
          fin_tforms = fopen(file_name_tforms, 'rt');
          transforms_tp.transforms = fscanf(fin_tforms, '%g', [6, inf]);
          fclose(fin_tforms);
          if(isempty(transforms_tp.transforms))
            transforms_tp.map_mask = [];
          else
            file_name_map = get_storage_file_name([file_name_suffix, '.map.tif']);
            if(exist(file_name_map, 'file')==2)
              fprintf('%s\n', file_name_map);
              transforms_tp.map_mask = imread(file_name_map);
            end
          end
        end
      else
        % Read in piecewise transformations as a MATLAB .mat file
        file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
          image_prefixes{plane_1, tile_a}, image_prefixes{plane_1, tile_b}, 'dt.');
        file_name = [file_name_suffix, '.mat'];
        try
          fprintf('%s\n', file_name);
          load2(file_name,'transforms_tp');
        catch %#ok<CTCH>
          transforms_tp = [];
        end
      end
      if(~isempty(transforms_tp) && ~isempty(transforms_tp.transforms) && ...
          ~isempty(transforms_tp.map_mask))
        matches = get_match_points_from_piecewise_affine(transforms_tp, ...
          fold_masks_1(tile_a).fold_mask, fold_masks_1(tile_b).fold_mask, 4);
      else
        matches.points_1 = [];
        matches.points_2 = [];
        matches.set_ids_1 = [];
        matches.set_ids_2 = [];
      end
      
      % Include Plan-B matches if asked to do so.
      if(isfield(dmesh_config, 'exceptions'))
        exception = ...
          get_exception(dmesh_config.exceptions, plane_1, 'plan_B_method');
      else
        exception = struct([]);
      end
      if(~(isfield(exception, 'plan_B_method') && ...
          strcmp(exception.plan_B_method, 'deformable_mesh')==0) && ...
          isfield(dmesh_config, 'plan_B_method') && ...
          strcmp(dmesh_config.plan_B_method, 'deformable_mesh')==1)
        file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
          image_prefixes{plane_1, tile_a}, image_prefixes{plane_1, tile_b}, 'dmcp.');
        file_name = [file_name_suffix, '.mat'];
        try
          fprintf('plan B: %s\n', file_name);
          plan_b = load2(file_name, 'match_points');
        catch
          plan_b.match_points.points_1 = [];
          plan_b.match_points.points_2 = [];
        end
        plan_b.match_points.set_ids_1 = [];
        plan_b.match_points.set_ids_2 = [];
        if(~isempty(plan_b.match_points.points_1) && ...
            ~isempty(plan_b.match_points.points_2))
          plan_b.match_points.points_1(1,:) = ...
            round(plan_b.match_points.points_1(1,:) + 1);
          plan_b.match_points.points_1(2,:) = ...
            round(plan_b.match_points.points_1(2,:) + 1);
          plan_b.match_points.points_2(1,:) = ...
            round(plan_b.match_points.points_2(1,:) + 1);
          plan_b.match_points.points_2(2,:) = ...
            round(plan_b.match_points.points_2(2,:) + 1);
          
          plan_b.match_points.set_ids_1 = fold_masks_1(tile_a).fold_mask(...
            sub2ind(size(fold_masks_1(tile_a).fold_mask), ...
            plan_b.match_points.points_1(2,:), ...
            plan_b.match_points.points_1(1,:)) );
          plan_b.match_points.set_ids_2 = fold_masks_1(tile_b).fold_mask(...
            sub2ind(size(fold_masks_1(tile_b).fold_mask), ...
            plan_b.match_points.points_2(2,:), ...
            plan_b.match_points.points_2(1,:)) );
        end
        matches.points_1 = [matches.points_1, plan_b.match_points.points_1];
        matches.points_2 = [matches.points_2, plan_b.match_points.points_2];
        matches.set_ids_1 = [matches.set_ids_1, plan_b.match_points.set_ids_1];
        matches.set_ids_2 = [matches.set_ids_2, plan_b.match_points.set_ids_2];
      end
      
      if(dmesh_config.is_verbose_figures)
        if(~isempty(matches.points_1))
          display_correspondence_points(images{plane_1, tile_a}, ...
            images{plane_1, tile_b}, matches.points_1, ...
            matches.points_2, -1);
          title(sprintf('plane: %d, tile 1: %d, tile 2: %d', ...
            case_ids_ref(plane_1), tile_a, tile_b));
        end
      end
      for patch_a = 1:n_patch(plane_1, tile_a)
        for patch_b = 1:n_patch(plane_1, tile_b)
          ids_to_add = matches.set_ids_1==patch_a & matches.set_ids_2==patch_b;
          match_points_a{patch_back_ref{plane_1,tile_a}(patch_a), ...
            patch_back_ref{plane_1,tile_b}(patch_b)} = [ ...
            match_points_a{patch_back_ref{plane_1,tile_a}(patch_a), ...
            patch_back_ref{plane_1,tile_b}(patch_b)}, ...
            matches.points_1(:, ids_to_add)];
          
          match_points_b{patch_back_ref{plane_1,tile_a}(patch_a), ...
            patch_back_ref{plane_1,tile_b}(patch_b)} = [ ...
            match_points_b{patch_back_ref{plane_1,tile_a}(patch_a), ...
            patch_back_ref{plane_1,tile_b}(patch_b)}, ...
            matches.points_2(:, ids_to_add)];
        end
      end
    end
  end
  
  %
  % deformable mesh transforms with other planes
  %
  print_debug('Correspondences across sections\n');
  for plane_2 = 1:size(tforms,1) % for each plane/section
    case_id_2 = case_ids_ref(plane_2);
    if(case_id_2<0 || plane_2==plane_1)
      continue;
    end
    if(nnz(adjacencies(plane_1,:,plane_2,:))<=0)
      continue;
    end
    if(dmesh_config.is_verbose)
      fprintf('plane1: %d, plane2: %d\n', case_id, case_id_2);
    end
    % load fold masks for plane2
    fold_masks_2 = struct('fold_mask', []);
    for tile = 1:size(tforms,2)
      if(isempty(tforms{plane_2,tile}))
        continue;
      end
      file_name_prefix =[fold_dir, image_prefixes{plane_2, tile}, '.fold_mask'];
      fold_masks_2(tile).fold_mask = load_fold_mask(file_name_prefix, config);
    end
    % generate correspondences and split according to folds
    print_debug('Generating correspondences ..\n');
    for tile_1 = 1:size(tforms,2)
      if(isempty(tforms{plane_1,tile_1}))
        continue;
      end
      for tile_2 = 1:size(tforms,2)
        if(isempty(tforms{plane_2,tile_2}))
          continue;
        end
        print_debug('tile pair: %d, %d\n', tile_1, tile_2);
        if(isfield(dmesh_config, 'input_tform_format') && ...
            strcmp(dmesh_config.input_tform_format, 'tif_txt'))
          % Read in piecewise transformations as a TIFF file of the map mask
          % and a text file of the affine transformations.
          file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
            image_prefixes{plane_1, tile_1}, image_prefixes{plane_2, tile_2}, 'dt.');
          file_name_tforms = get_storage_file_name([file_name_suffix, '.tforms.txt']);
          if(exist(file_name_tforms, 'file')~=2)
            continue;
          end
          transforms_tp = [];
          fprintf('%s\n', file_name_tforms);
          fin_tforms = fopen(file_name_tforms, 'rt');
          transforms_tp.transforms = fscanf(fin_tforms, '%g', [6, inf]);
          fclose(fin_tforms);
          if(isempty(transforms_tp.transforms))
            transforms_tp.map_mask = [];
            continue;
          end
          file_name_map = get_storage_file_name([file_name_suffix, '.map.tif']);
          if(exist(file_name_map, 'file')~=2)
            continue;
          end
          fprintf('%s\n', file_name_map);
          transforms_tp.map_mask = imread(file_name_map);
        else
          file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
            image_prefixes{plane_1, tile_1}, image_prefixes{plane_2, tile_2}, 'dt.');
          file_name = [file_name_suffix, '.mat'];
          try
            fprintf('%s\n', file_name);
            load2(file_name,'transforms_tp');
          catch
            continue;
          end
        end
        if(isempty(transforms_tp.map_mask))
          continue;
        end
        matches = get_match_points_from_piecewise_affine(...
          transforms_tp, fold_masks_1(tile_1).fold_mask, ...
          fold_masks_2(tile_2).fold_mask, 4);
        
        % include exception for minimum number of matches
        if(isfield(dmesh_config, 'exceptions'))
          exception = ...
            get_exception(dmesh_config.exceptions, case_id, 'min_n_matches');
          if(isfield(exception, 'min_n_matches') && ...
              ~isempty(exception.min_n_matches))
            [junk, n_m1] = count_row_occurence(matches.set_ids_1');
            [junk, n_m2] = count_row_occurence(matches.set_ids_2');
            n_m = min([n_m1; n_m2]);
            if(n_m<=dmesh_config.min_n_matches && ...
                n_m>exception.min_n_matches)
              matches.points_1 = [matches.points_1, matches.points_1];
              matches.points_2 = [matches.points_2, matches.points_2];
              matches.set_ids_1 = [matches.set_ids_1, matches.set_ids_1];
              matches.set_ids_2 = [matches.set_ids_2, matches.set_ids_2];
            end
          end
        end
        
        % Include Plan-B matches if asked to do so.
        if(isfield(dmesh_config, 'exceptions'))
          exception = ...
            get_exception(dmesh_config.exceptions, case_id, 'plan_B_method');
        else
          exception = struct([]);
        end
        if(~(isfield(exception, 'plan_B_method') && ...
            strcmp(exception.plan_B_method, 'deformable_mesh')==0) && ...
            isfield(dmesh_config, 'plan_B_method') && ...
            strcmp(dmesh_config.plan_B_method, 'deformable_mesh')==1)
          file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
            image_prefixes{plane_1, tile_1}, image_prefixes{plane_2, tile_2}, 'dmcp.');
          file_name = [file_name_suffix, '.mat'];
          try
            fprintf('plan B: %s\n', file_name);
            plan_b = load2(file_name, 'match_points');
          catch
            plan_b.match_points.points_1 = [];
            plan_b.match_points.points_2 = [];
          end
          plan_b.match_points.set_ids_1 = [];
          plan_b.match_points.set_ids_2 = [];
          if(~isempty(plan_b.match_points.points_1) && ...
              ~isempty(plan_b.match_points.points_2))
            plan_b.match_points.points_1(1,:) = ...
              round(plan_b.match_points.points_1(1,:) + 1);
            plan_b.match_points.points_1(2,:) = ...
              round(plan_b.match_points.points_1(2,:) + 1);
            plan_b.match_points.points_2(1,:) = ...
              round(plan_b.match_points.points_2(1,:) + 1);
            plan_b.match_points.points_2(2,:) = ...
              round(plan_b.match_points.points_2(2,:) + 1);
            
            plan_b.match_points.set_ids_1 = fold_masks_1(tile_1).fold_mask(...
              sub2ind(size(fold_masks_1(tile_1).fold_mask), ...
              plan_b.match_points.points_1(2,:), ...
              plan_b.match_points.points_1(1,:)) );
            plan_b.match_points.set_ids_2 = fold_masks_1(tile_2).fold_mask(...
              sub2ind(size(fold_masks_1(tile_2).fold_mask), ...
              plan_b.match_points.points_2(2,:), ...
              plan_b.match_points.points_2(1,:)) );
          end
          matches.points_1 = [matches.points_1, plan_b.match_points.points_1];
          matches.points_2 = [matches.points_2, plan_b.match_points.points_2];
          matches.set_ids_1 = [matches.set_ids_1, plan_b.match_points.set_ids_1];
          matches.set_ids_2 = [matches.set_ids_2, plan_b.match_points.set_ids_2];
        end
        
        if(dmesh_config.is_verbose_figures)
          if(~isempty(matches.points_1))
            display_correspondence_points(images{plane_1, tile_1}, ...
              images{plane_2, tile_2}, matches.points_1, ...
              matches.points_2, -1);
            title(sprintf('plane 1: %d, plane 2: %d, tile 1: %d, tile 2: %d', ...
              case_ids_ref(plane_1), case_ids_ref(plane_2), tile_1, tile_2));
          end
        end
        for patch_a = 1:n_patch(plane_1, tile_1)
          for patch_b = 1:n_patch(plane_2, tile_2)
            ids_to_add = matches.set_ids_1==patch_a & matches.set_ids_2==patch_b;
            match_points_a{patch_back_ref{plane_1,tile_1}(patch_a), ...
              patch_back_ref{plane_2,tile_2}(patch_b)} = [ ...
              match_points_a{patch_back_ref{plane_1,tile_1}(patch_a), ...
              patch_back_ref{plane_2,tile_2}(patch_b)}, ...
              matches.points_1(:, ids_to_add)];
            
            match_points_b{patch_back_ref{plane_1,tile_1}(patch_a), ...
              patch_back_ref{plane_2,tile_2}(patch_b)} = [ ...
              match_points_b{patch_back_ref{plane_1,tile_1}(patch_a), ...
              patch_back_ref{plane_2,tile_2}(patch_b)}, ...
              matches.points_2(:, ids_to_add)];
          end
        end
      end
    end
  end
  print_debug('---\n');
end

check_for_dir('~/temp/LM_exp/');
save(['~/temp/LM_exp/global_align_match_points.', ...
  num2str(stack_config.case_ids(1)), '.', ...
  num2str(stack_config.case_ids(end)), '.mat'], ...
  'match_points_a', 'match_points_b');

load(['~/temp/LM_exp/global_align_match_points.', ...
  num2str(stack_config.case_ids(1)), '.', ...
  num2str(stack_config.case_ids(end)), '.mat'], ...
  'match_points_a', 'match_points_b');

% make matches symmetric
for patch_id_1 = 1:n_patch_global
  for patch_id_2 = 1:n_patch_global
    if(size(match_points_a{patch_id_1, patch_id_2},2)<=dmesh_config.min_n_matches)
      match_points_a{patch_id_1, patch_id_2} = [];
      match_points_b{patch_id_1, patch_id_2} = [];
    end
  end
end

for patch_id_1 = 1:n_patch_global
  for patch_id_2 = 1:patch_id_1
    if(size(match_points_a{patch_id_1, patch_id_2},2)<=dmesh_config.min_n_matches && ...
        size(match_points_b{patch_id_2, patch_id_1},2)<=dmesh_config.min_n_matches )
      m_p_a = [];
      m_p_b = [];
    else
      m_p_a = [match_points_a{patch_id_1, patch_id_2}, ...
        match_points_b{patch_id_2, patch_id_1}];
      m_p_b = [match_points_b{patch_id_1, patch_id_2}, ...
        match_points_a{patch_id_2, patch_id_1}];
    end
    
    match_points_a{patch_id_1, patch_id_2} = m_p_a;
    match_points_b{patch_id_2, patch_id_1} = m_p_a;
    
    match_points_b{patch_id_1, patch_id_2} = m_p_b;
    match_points_a{patch_id_2, patch_id_1} = m_p_b;
    
    %     print_debug('match_points_a{%d, %d}:\n', patch_id_1, patch_id_2);
    %     disp_debug(match_points_a{patch_id_1, patch_id_2});
    %     print_debug('match_points_b{%d, %d}:\n', patch_id_1, patch_id_2);
    %     disp_debug(match_points_b{patch_id_1, patch_id_2});
    %     print_debug('match_points_a{%d, %d}:\n', patch_id_2, patch_id_1);
    %     disp_debug(match_points_a{patch_id_2, patch_id_1});
    %     print_debug('match_points_b{%d, %d}:\n', patch_id_2, patch_id_1);
    %     disp_debug(match_points_b{patch_id_2, patch_id_1});
  end
end

% build a graph with number of matches as weights and find the largest
% connected component
adjacency_matrix = sparse(zeros(n_patch_global));
for patch_id_1 = 1:n_patch_global
  for patch_id_2 = 1:n_patch_global
    if(~isempty(match_points_a{patch_id_1, patch_id_2}))
      adjacency_matrix(patch_id_1, patch_id_2) = ...
        size(match_points_a{patch_id_1, patch_id_2},2);
    end
  end
end
print_debug('adjacency_matrix:\n');
disp_debug(adjacency_matrix);
connected_component_partition = double(get_connected_components(adjacency_matrix, 3));
print_debug('connected_component_partition:\n');
disp_debug(connected_component_partition);
largest_component_id = mode(connected_component_partition);
print_debug('largest_component_id:\n');
disp_debug(largest_component_id);
is_member_largest_component = connected_component_partition==largest_component_id;
print_debug('is_member_largest_component:\n');
disp_debug(is_member_largest_component);
n_member_largest_component = sum(is_member_largest_component);
print_debug('n_member_largest_component:\n');
disp_debug(n_member_largest_component);

% truncate the matches to those patches belonging to the largest component
match_points_1_a{n_member_largest_component, n_member_largest_component} = [];
match_points_1_b{n_member_largest_component, n_member_largest_component} = [];
member_id_1 = 0;
new_to_original_map = zeros(1, n_member_largest_component);
original_to_new_map = zeros(1, n_patch_global);
for patch_id_1 = 1:n_patch_global
  if(~is_member_largest_component(patch_id_1))
    continue;
  end
  member_id_1 = member_id_1 + 1;
  new_to_original_map(member_id_1) = patch_id_1;
  original_to_new_map(patch_id_1) = member_id_1;
  for patch_id_2 = 1:patch_id_1
    if(~is_member_largest_component(patch_id_2))
      continue;
    end
    match_points_1_a{original_to_new_map(patch_id_1), ...
      original_to_new_map(patch_id_2)} = match_points_a{patch_id_1, patch_id_2};
    match_points_1_b{original_to_new_map(patch_id_1), ...
      original_to_new_map(patch_id_2)} = match_points_b{patch_id_1, patch_id_2};
    
    match_points_1_a{original_to_new_map(patch_id_2), ...
      original_to_new_map(patch_id_1)} = match_points_a{patch_id_2, patch_id_1};
    match_points_1_b{original_to_new_map(patch_id_2), ...
      original_to_new_map(patch_id_1)} = match_points_b{patch_id_2, patch_id_1};
  end
end

% get the patch with maximum number of matches and assign identity
% transformation to it
max_n_match = 0;
max_n_match_patch_id = -1;
for patch_id_1 = 1:size(match_points_1_a, 1)
  n_match = 0;
  for patch_id_2 = 1:size(match_points_1_a, 2)
    n_match = n_match + size(match_points_1_a{patch_id_1, patch_id_2},2);
  end
  if(n_match>max_n_match)
    max_n_match = n_match;
    max_n_match_patch_id = patch_id_1;
  end
end
if(dmesh_config.is_verbose)
  fprintf('Assigning identity transform to patch %d\n', max_n_match_patch_id);
end

% compute affine transformations for the patches belonging to largest
% component
alpha = ones(size(match_points_1_a));

if(~isfield(dmesh_config, 'optimization'))
  dmesh_config.optimization.method = 'least_squares';
end
switch(dmesh_config.optimization.method)
  case 'least_squares'
    if(dmesh_config.is_verbose)
      fprintf('Performing Least-Squares based optimization ...');
    end
    [a,b] = getSIFTmodel(match_points_1_a, match_points_1_b, ...
      alpha, [], max_n_match_patch_id);
    if(dmesh_config.is_verbose)
      fprintf('done\n');
    end
  case 'levenberg_marquadt_rigid'
    if(dmesh_config.is_verbose)
      fprintf('Performing Levenberg-Marquadt based optimization ...');
    end
    dmesh_config.optimization.levenberg_marquadt_rigid.is_verbose = ...
      dmesh_config.is_verbose;
    dmesh_config.optimization.levenberg_marquadt_rigid.is_verbose_figures = ...
      dmesh_config.is_verbose_figures;
    [a,b] = ...
      get_affine_levenberg_marquadt_rigid(...
      match_points_1_a, match_points_1_b, max_n_match_patch_id, ...
      dmesh_config.optimization.levenberg_marquadt_rigid);
    if(dmesh_config.is_verbose)
      fprintf('done\n');
    end
  case 'least_squares_rigidified'
    if(dmesh_config.is_verbose)
      fprintf('Performing Least-Squares alogn with rigid based optimization ...');
    end
    [a,b] = least_squares_rigidified(match_points_1_a, match_points_1_b, ...
      alpha, max_n_match_patch_id);
    if(dmesh_config.is_verbose)
      fprintf('done\n');
    end
end
D = [a, b];

if(dmesh_config.is_verbose)
  fprintf('Collecting together the global transformations ...');
end
for plane=1:size(tforms,1)
  for tile=1:size(tforms,2)
    cnt = patch_back_ref{plane,tile};
    if(isempty(cnt))
      continue;
    end
    transforms(plane,tile).transform = [];
    for c = cnt
      print_debug('plane %d, tile %d, patch %d\n', case_ids_ref(plane), tile, c);
      if(~is_member_largest_component(c))
        % this patch did not have adequate matches so couldn't compute
        % transform. Copy transform from some other patch in the image
        i = find(is_member_largest_component(patch_back_ref{plane,tile}));
        if(~isempty(i))
          if(dmesh_config.is_verbose)
            fprintf('Found a patch that did not have its own transform. Therefore, \n');
            fprintf('copying a transform from another patch in this tile.\n');
          end
          copy_patch_id_1 = original_to_new_map(patch_back_ref{plane,tile}(i(1)));
          transforms(plane,tile).transform(:,end+1) = ...
            [reshape(D{copy_patch_id_1}, [1 4]), ...
            D{copy_patch_id_1+n_member_largest_component}']';
        else
          if(dmesh_config.is_verbose)
            fprintf('Found a patch that did not have its own transform. Could not\n');
            fprintf('even find any other patch in this tile with a valid transform.\n');
            fprintf('Simply storing identity transform.\n');
          end
          transforms(plane,tile).transform(:,end+1) = [1 0 0 1 0 0];
        end
      else
        c_1 = original_to_new_map(c);
        transforms(plane,tile).transform(:,end+1) = ...
          [reshape(D{c_1}, [1 4]), D{c_1+n_member_largest_component}']';
      end
    end
  end
end
if(dmesh_config.is_verbose)
  fprintf('done\n');
end

%Save files
if(dmesh_config.is_verbose)
  fprintf('Saving transformations ... ');
end
tform_dir = [get_region_dir(config), config.deformable_mesh.dir];
if(~isempty(transforms))
  for plane=1:size(transforms,1)
    for tile=1:size(transforms,2)
      check_for_dir([tform_dir, image_sub_dirs{plane, tile}]);
      file_name = [tform_dir, image_prefixes{plane, tile}, ...
        '.patchwise_affine_transforms.mat'];
      transform=transforms(plane, tile); %#ok<NASGU>
      save2(file_name,'transform');
    end
  end
end
if(dmesh_config.is_verbose)
  fprintf('done\n');
end

return
end
