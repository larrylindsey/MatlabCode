function align_stack_dmesh_affine_in_plane(config)
%  align_stack_dmesh_affine_tile_patches(config)
% Compute affine transformations based on correspondences from deformable
% mesh.
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
dmesh_config = config.align.in_section_align.deformable_mesh;

dmesh_dir = get_deformable_mesh_dir(config);

if(~isfield(dmesh_config, 'is_verbose'))
  dmesh_config.is_verbose = true;
end
if(~isfield(dmesh_config, 'is_verbose_figures'))
  dmesh_config.is_verbose_figures = false;
end
if(dmesh_config.is_verbose)
  fprintf('START: align_stack_dmesh_affine_in_plane\n');
end

%%%%%%%% Load correspondences from deformable mesh alignments %%%%%%
% An affine transformation is computed for patch in the global reference
% frame.
% Each patch has a within-tile and global id. This mapping is maintained in
% patch_back_ref{tile} =
%     [global_id_of_patch_1, global_id_of_patch_2, ... ]
if(dmesh_config.is_verbose)
  fprintf('Initializing structures .. \n');
end
for plane = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(plane);
  if(dmesh_config.is_verbose)
    fprintf('case_id: %d\n', case_id);
  end
  
  [image_prefixes, image_sub_dirs, is_to_be_processed] = ...
    get_image_prefixes_subdirs(config, case_id);

  if(max(is_to_be_processed)==0)
    if(dmesh_config.is_verbose)
      fprintf('All tiles in this section have is_to_be_processed=false, skipping this section\n');
    end
    continue;
  end
  
  if(dmesh_config.is_verbose_figures)
    images = get_image_from_stack(config, case_id);
  end
  
  n_patch_global = 0;
  patch_back_ref = cell(1, length(image_prefixes));
  n_patch = zeros(1, length(image_prefixes));
  
  for tile = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile}))
      continue;
    end
    n_patch(tile) = 1;
    patch_back_ref{tile} = ...
      n_patch_global+1 : n_patch_global+n_patch(tile);
    n_patch_global = n_patch_global + n_patch(tile);
  end
  print_debug('n_patch\n');
  disp_debug(n_patch);
  if(dmesh_config.is_verbose)
    fprintf('done\n');
  end
  
  match_points_a = cell(n_patch_global, n_patch_global);
  match_points_b = cell(n_patch_global, n_patch_global);
  
  %
  % deformable mesh transforms for within this plane
  %
  % generate correspondences and split according to folds
  print_debug('Correspondences within section\n');
  for tile_a = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_a}))
      continue;
    end
    for tile_b = tile_a+1:length(image_prefixes)
      if(isempty(image_prefixes{tile_b}))
        continue;
      end
      print_debug('tile pair: %d, %d\n', tile_a, tile_b);
      print_debug('%s\n%s\n', image_prefixes{tile_a}, ...
        image_prefixes{tile_b});
      transforms_tp = [];
      if(isfield(dmesh_config, 'input_tform_format') && ...
          strcmp(dmesh_config.input_tform_format, 'tif_txt'))
        fprintf('Trying to read dmesh .tif and .txt files\n');
        % Read in piecewise transformations as a TIFF file of the map mask
        % and a text file of the affine transformations.
        file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
          image_prefixes{tile_a}, image_prefixes{tile_b}, 'dt.');
        file_name_tforms = get_storage_file_name([file_name_suffix, '.tforms.txt']);
        if(exist(file_name_tforms, 'file')==2)
          fprintf('transform .txt exists\n');
          transforms_tp = [];
          fin_tforms = fopen(file_name_tforms, 'rt');
          transforms_tp.transforms = fscanf(fin_tforms, '%g', [6, inf]);
          fclose(fin_tforms);
          if(isempty(transforms_tp.transforms))
            fprintf('transform .txt is empty\n');
            transforms_tp.map_mask = [];
            fprintf('transform file txt:\n%s\n', ...
              get_storage_file_name(file_name_tforms));
          else
            fprintf('Successfully read transform .txt file. Trying to read .tif file\n');
            file_name_map = get_storage_file_name([file_name_suffix, '.map.tif']);
            if(exist(file_name_map, 'file')==2)
              fprintf('Transform .tif file exists\n');
              transforms_tp.map_mask = imread(file_name_map);
              fprintf('transform file tif:\n%s\n', ...
                get_storage_file_name(file_name_map));
            else
              fprintf('Transform .tif file does not exist\n');
            end
          end
        else
          fprintf('transform .txt file does not exist\n');
        end
      else
        fprintf('Trying to read dmesh .mat file\n');
        % Read in piecewise transformations as a MATLAB .mat file
        file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
          image_prefixes{tile_a}, image_prefixes{tile_b}, 'dt.');
        file_name = [file_name_suffix, '.mat'];
        try
          load2(file_name,'transforms_tp');
          fprintf('transform file mat:\n%s\n', ...
            get_storage_file_name(file_name));
          fprintf('Successfully read dmesh .mat file\n');
        catch %#ok<CTCH>
          fprintf('Could not read dmesh .mat file\n%s\n', ...
            get_storage_file_name(file_name));
          transforms_tp = [];
        end
      end
      if(~isempty(transforms_tp) && ~isempty(transforms_tp.transforms) && ...
          ~isempty(transforms_tp.map_mask))
        matches = get_match_points_from_piecewise_affine(transforms_tp);
      else
        matches.points_1 = [];
        matches.points_2 = [];
      end
      
      % Include Plan-B matches if asked to do so.
      if(isfield(dmesh_config, 'plan_B_method') && ...
          strcmp(dmesh_config.plan_B_method, 'deformable_mesh')==1)
        file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
          image_prefixes{tile_a}, image_prefixes{tile_b}, 'dmcp.');
        file_name = [file_name_suffix, '.mat'];
        try
          plan_b = load2(file_name, 'match_points');
          fprintf('transform file mat:\n%s\n', ...
            get_storage_file_name(file_name));
        catch %#ok<CTCH>
          plan_b.match_points.points_1 = [];
          plan_b.match_points.points_2 = [];
        end
        
        matches.points_1 = [matches.points_1, plan_b.match_points.points_1];
        matches.points_2 = [matches.points_2, plan_b.match_points.points_2];
      end
      
      if(dmesh_config.is_verbose_figures)
        if(~isempty(matches.points_1))
          display_correspondence_points(images{plane, tile_a}, ...
            images{plane, tile_b}, matches.points_1, ...
            matches.points_2, -1);
          title(sprintf('plane: %d, tile 1: %d, tile 2: %d', case_id, tile_a, tile_b));
        end
      end
      for patch_a = 1:n_patch( tile_a)
        for patch_b = 1:n_patch( tile_b)
          match_points_a{patch_back_ref{tile_a}(patch_a), ...
            patch_back_ref{tile_b}(patch_b)} = [ ...
            match_points_a{patch_back_ref{tile_a}(patch_a), ...
            patch_back_ref{tile_b}(patch_b)}, ...
            matches.points_1];
          
          match_points_b{patch_back_ref{tile_a}(patch_a), ...
            patch_back_ref{tile_b}(patch_b)} = [ ...
            match_points_b{patch_back_ref{tile_a}(patch_a), ...
            patch_back_ref{tile_b}(patch_b)}, ...
            matches.points_2];
        end
      end
    end
  end
  
  
  % make matches symmetric
  for patch_id_1 = 1:n_patch_global
    for patch_id_2 = 1:patch_id_1
      if(size(match_points_a{patch_id_1, patch_id_2},2)<1 && ...
          size(match_points_b{patch_id_2, patch_id_1},2)<1 )
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
      
      print_debug('match_points_a{%d, %d}:\n', patch_id_1, patch_id_2);
      disp_debug(match_points_a{patch_id_1, patch_id_2});
      print_debug('match_points_b{%d, %d}:\n', patch_id_1, patch_id_2);
      disp_debug(match_points_b{patch_id_1, patch_id_2});
      print_debug('match_points_a{%d, %d}:\n', patch_id_2, patch_id_1);
      disp_debug(match_points_a{patch_id_2, patch_id_1});
      print_debug('match_points_b{%d, %d}:\n', patch_id_2, patch_id_1);
      disp_debug(match_points_b{patch_id_2, patch_id_1});
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
  connected_component_partition = double(get_connected_components(adjacency_matrix, 1));
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
  match_points_1_a = ...
    cell(n_member_largest_component, n_member_largest_component);
  match_points_1_b = ...
    cell(n_member_largest_component, n_member_largest_component);
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
  
  % compute affine transformations for the patches belonging to largest
  % component
  alpha = zeros(size(match_points_1_a));
  for patch_id_1 = 1:size(match_points_1_a,1)
    for patch_id_2 = 1:size(match_points_1_a,2)
      if(size(match_points_1_a, 2)>0)
        alpha(patch_id_1, patch_id_2) = 1;
      end
    end
  end
  
  D = compute_affine_transform(dmesh_config, match_points_1_a, ...
    match_points_1_b, alpha);
  is_valid_transform = get_transform_validity(D);
  if(~is_valid_transform)
    % There are some false correspondences - try to remove them
    for i = 1:numel(alpha)
      fprintf(['Found suspicious transform.', ...
        'Trying to remove false correspondence.\n']);
      if(alpha(i)==0)
        continue;
      end
      i1 = mod(i, size(alpha,1))+1;
      i2 = ceil(i/size(alpha,1));
      alpha(i1,i2) = 0;
      alpha(i2,i1) = 0;
      
      D = compute_affine_transform(dmesh_config, match_points_1_a, ...
        match_points_1_b, alpha);
      is_valid_transform = get_transform_validity(D);
      if(is_valid_transform)
        break;
      end
      alpha(i1,i2) = 1;
      alpha(i2,i1) = 1;
    end
  end
  if(~is_valid_transform)
    warning('ALIGN_STACK_DMESH:INVALID_TFORM', 'Suspicious transformation');
  end
  
  if(dmesh_config.is_verbose)
    fprintf('Collecting together the global transformations ...\n');
  end
  transforms = [];
  for tile=1:length(image_prefixes)
    cnt = patch_back_ref{tile};
    if(isempty(cnt))
      continue;
    end
    transforms(tile).transform = []; %#ok<AGROW>
    for c = cnt
      print_debug('plane %d, tile %d, patch %d\n', case_id, tile, c);
      print_debug('%s\n', image_prefixes{tile});
      if(~is_member_largest_component(c))
        % this patch did not have adequate matches so couldn't compute
        % transform. Copy transform from some other patch in the image
        i = find(is_member_largest_component(patch_back_ref{tile}));
        if(~isempty(i))
          if(dmesh_config.is_verbose)
            fprintf('Found a patch that did not have its own transform. Therefore, \n');
            fprintf('copying a transform from another patch in this tile.\n');
          end
          copy_patch_id_1 = original_to_new_map(patch_back_ref{tile}(i(1)));
          transforms(tile).transform(:,end+1) = ...
            [reshape(D{copy_patch_id_1}, [1 4]), ...
            D{copy_patch_id_1+n_member_largest_component}']'; %#ok<AGROW>
        else
          if(dmesh_config.is_verbose)
            fprintf('Found a patch that did not have its own transform. Could not\n');
            fprintf('even find any other patch in this tile with a valid transform.\n');
            fprintf('Simply storing identity transform.\n');
          end
          transforms(tile).transform(:,end+1) = [1 0 0 1 0 0]; %#ok<AGROW>
        end
      else
        c_1 = original_to_new_map(c);
        transforms(tile).transform(:,end+1) = ...
          [reshape(D{c_1}, [1 4]), D{c_1+n_member_largest_component}']'; %#ok<AGROW>
      end
    end
  end
  if(dmesh_config.is_verbose)
    fprintf('done\n');
  end
  
  %
  % Save files
  %
  if(dmesh_config.is_verbose)
    fprintf('Saving transformations ...\n');
  end
  % Generate a substring for the file name that accounts for the images
  % involved in the joint alignment
  image_set_string = '.';
  for tile = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile}))
      continue;
    end
    image_set_string = [image_set_string, image_prefixes{tile}]; %#ok<AGROW>
  end
  image_set_string = strrep(image_set_string, '/', '_');
  image_set_string = strrep(image_set_string, '\', '_');
  if(~isempty(transforms))
    for tile=1:size(transforms,2)
      check_for_dir([dmesh_dir, image_sub_dirs{tile}]);
      file_name = [dmesh_dir, image_prefixes{tile}, image_set_string, ...
        '.dmesh_affine.in_plane.mat'];
      transform_p = transforms(tile).transform; %#ok<NASGU>
      save2(file_name,'transform_p');
    end
  else
    warning('transforms found to be empty'); %#ok<WNTAG>
  end
  if(dmesh_config.is_verbose)
    fprintf('done\n');
  end
  print_debug('---\n');
end

return
end

function is_valid_transform = get_transform_validity(D)
is_valid_transform = true;
for t = 1:length(D)/2
  if(sum(isnan(D{t}(:)))>0)
    is_valid_transform = false;
    fprintf('Transform has nan\n');
    break;
  end
  if(sum(abs(D{t}(:))>2)>0)
    is_valid_transform = false;
    fprintf('Transform has abs value greater than 2\n');
    display(D{t});
    break;
  end
  det_t = det(D{t});
  if(det_t > 1.2 || det_t < 0.8)
    is_valid_transform = false;
    fprintf('Suspicious determinant: %d\n', det_t);
    break;
  end
  [e,l] = eig(D{t});
  rel_eig_val_diff = 2*abs(l(1,1)-l(2,2))/(abs(l(1,1))+abs(l(2,2)));
  if(rel_eig_val_diff >0.5)
    is_valid_transform = false;
    fprintf('Eigen values are quite different: %d\n', rel_eig_val_diff);
    break;
  end
end
return
end

function D = compute_affine_transform(dmesh_config, match_points_1_a, ...
  match_points_1_b, alpha)
if(~isfield(dmesh_config, 'optimization'))
  dmesh_config.optimization.method = 'least_squares';
end
switch(dmesh_config.optimization.method)
  case 'least_squares'
    if(dmesh_config.is_verbose)
      fprintf('Performing Least-Squares based optimization ...');
    end
    [junk, fix_id] = max(sum(alpha, 2));
    [a,b] = getSIFTmodel(match_points_1_a, match_points_1_b, alpha, [], fix_id);
    if(dmesh_config.is_verbose)
      fprintf('done\n');
    end
  case 'levenberg_marquadt_rigid'
    if(dmesh_config.is_verbose)
      fprintf('Performing Levenberg-Marquadt based optimization ...');
    end
    [a,b] = ...
      get_affine_levenberg_marquadt_rigid(match_points_1_a, match_points_1_b);
    if(dmesh_config.is_verbose)
      fprintf('done\n');
    end
  otherwise
    error('Optimization method not understood.\n');
end
D = [a, b];
end
