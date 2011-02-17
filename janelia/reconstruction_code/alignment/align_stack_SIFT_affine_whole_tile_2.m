function align_stack_SIFT_affine_whole_tile_2(config)
%  align_stack_dmesh_affine_tile_patches(config)
% Compute affine transformations based on correspondences from SIFT feature
% points. Pairs of images are aligned apriori. The obtained inlying
% correspondences are used to compute global joint alignment of all images.
%
% Input:
%   config    config datastructure of the reconstruction
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  010709  init. code
% v1  030209  code modified from align_stack_dmesh_affine_tile_pieces.m
%

stack_config = config.stack;
sift_config = config.align.global_align.SIFT;

image_dir = get_stack_dir(config);
sift_dir = get_sift_dir(config);

if( ~isfield(stack_config, 'image_structure') || isempty(stack_config.image_structure) )
  error('xml file from TrakEM not specified. Could not load tile images and initial transformations.');
end;

if(~isfield(sift_config, 'is_verbose'))
  sift_config.is_verbose = true;
end
if(~isfield(sift_config, 'is_verbose_figures'))
  sift_config.is_verbose_figures = false;
end
if(sift_config.is_verbose)
  fprintf('Global alignment based on SIFT correspondences\n');
end

s = xmlread([image_dir, stack_config.image_structure]);
xmlwrite('~tmp.xml', s);
xmlstr = fileread('~tmp.xml');
s = xml_parseany(xmlstr);
s_t2_layer = s.t2_layer_set{1}.t2_layer;

if(sift_config.is_verbose)
  fprintf('Collecting stack images and initial trakEM transformations ..');
end
if(sift_config.is_verbose_figures)
  images = {};
end
image_prefixes = {};
image_sub_dirs = {};
tforms={};
sizes={};
n_patch_global = 0;
patch_back_ref = {};
case_ids_ref = -ones(1,length(s_t2_layer)); % record of case_id of each plane to be used when saving
for i = 1:length(s_t2_layer)
  case_id=str2double(s_t2_layer{i}.ATTRIBUTE.z);
  if(ismember(case_id,stack_config.case_ids) || isempty(stack_config.case_ids))
    case_ids_ref(i) = case_id;
    for j = 1:length(s_t2_layer{i}.t2_patch)
      % remove the directory and extension from the file_path and store
      % in image_prefix
      dir_length = strfind(s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path, '/');
      dir_length = [zeros(1, 3-length(dir_length)), dir_length];

      % File path in .xml file is
      % <stack_directory>/<grid_id>/<section_id>/<image_file_name>
      if(sift_config.is_verbose_figures)
        images{i,j} = im2double(imread(...
          [image_dir, s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path(...
          dir_length(end-2)+1:end)]));
        if(size(images{i,j},3)>1)
          images{i,j} = rgb2gray(images{i,j});
        end
        if(max(images{i,j}(:))>1)
          images{i,j} = images{i,j}/max(images{i,j}(:));
        end
      end
      image_prefixes{i,j} = ...
        s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path(dir_length(end-2)+1:end-4);
      image_sub_dirs{i,j} = ...
        s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path(dir_length(end-2)+1:dir_length(end));

      n_patch_global = n_patch_global + 1;
      patch_back_ref{i,j} = n_patch_global;

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
if(sift_config.is_verbose)
  fprintf('done.\n');
end


%Generate adjacency matrix from trakEM2 data
if(sift_config.trakEM2.getAdjacency && ~isempty(tforms))
  if(sift_config.is_verbose)
    fprintf('Computing tile adjacency matrix ..\n');
  end
  adjacencies=getSIFTalpha(tforms,sizes);
  if(sift_config.is_verbose)
    fprintf('done.\n');
  end
else
  adjacencies=[];
end

if(sift_config.save)
  if(~isempty(adjacencies))
    file_name=[sift_dir, 'adjacencies.mat'];
    save2(file_name, 'adjacencies');
  end
end

%%%%%%%% Load correspondences from SIFT alignments %%%%%%
match_points_a{n_patch_global, n_patch_global} = [];
match_points_b{n_patch_global, n_patch_global} = [];
match_error = zeros(n_patch_global);

for plane_1 = 1:size(tforms,1) % for each plane/section
  case_id = case_ids_ref(plane_1);
  if(case_id<0)
    continue;
  end
  if(sift_config.is_verbose)
    fprintf('plane1: %d\n', case_id);
  end

  %
  % SIFT inlier matches for within this plane
  %
  % generate correspondences and split according to folds
  print_debug('Correspondences within section\n');
  for tile_a = 1:size(tforms,2)
    if(isempty(tforms{plane_1,tile_a}))
      continue;
    end
    for tile_b = 1:size(tforms,2)
      if(isempty(tforms{plane_1,tile_b}))
        continue;
      end
      file_name_suffix = get_file_name_from_tuple(sift_dir, ...
        image_prefixes{plane_1, tile_a}, image_prefixes{plane_1, tile_b}, 'st.');
      file_name = [file_name_suffix, '.mat'];
      try
        load2(file_name,'transforms_tp', 'err', 'match_points');
      catch
        continue;
      end
      if(isempty(match_points(1,2).points_1))
        continue;
      end
      print_debug('tile pair: %d, %d\n', tile_a, tile_b);
      if(sift_config.is_verbose_figures)
        display_correspondence_points(images{plane_1, tile_a}, ...
          images{plane_1, tile_b}, match_points(1,2).points_1, ...
          match_points(1,2).points_2, -1);
        title(sprintf('plane: %d, tile 1: %d, tile 2: %d', plane_1, tile_a, tile_b));
      end
      match_points_a{patch_back_ref{plane_1,tile_a}, ...
        patch_back_ref{plane_1,tile_b}} = [ ...
        match_points_a{patch_back_ref{plane_1,tile_a}, ...
        patch_back_ref{plane_1,tile_b}}, match_points(1,2).points_1];

      match_points_b{patch_back_ref{plane_1,tile_a}, ...
        patch_back_ref{plane_1,tile_b}} = [ ...
        match_points_b{patch_back_ref{plane_1,tile_a}, ...
        patch_back_ref{plane_1,tile_b}}, match_points(1,2).points_2];
      
      match_error(patch_back_ref{plane_1,tile_a}, patch_back_ref{plane_1,tile_b}) = ...
        err;
    end
  end

  %
  % SIFT correspondences with other planes
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
    if(sift_config.is_verbose)
      fprintf('plane1: %d, plane2: %d\n', case_id, case_id_2);
    end

    % correspondences
    print_debug('Collecting correspondences ..\n');
    for tile_1 = 1:size(tforms,2)
      if(isempty(tforms{plane_1,tile_1}))
        continue;
      end
      for tile_2 = 1:size(tforms,2)
        if(isempty(tforms{plane_2,tile_2}))
          continue;
        end
        file_name_suffix = get_file_name_from_tuple(sift_dir, ...
          image_prefixes{plane_1, tile_1}, image_prefixes{plane_2, tile_2}, 'st.');
        file_name = [file_name_suffix, '.mat'];
        try
          load2(file_name,'transforms_tp', 'err', 'match_points');
        catch
          continue;
        end
        if(isempty(match_points(1,2).points_1))
          continue;
        end
        print_debug('tile pair: %d, %d\n', tile_1, tile_2);
        if(sift_config.is_verbose_figures)
          display_correspondence_points(images{plane_1, tile_1}, ...
            images{plane_2, tile_2}, match_points(1,2).points_1, ...
            match_points(1,2).points_2, -1);
          title(sprintf('plane 1: %d, plane 2: %d, tile 1: %d, tile 2: %d', ...
            plane_1, plane_2, tile_1, tile_2));
        end
        match_points_a{patch_back_ref{plane_1,tile_1}, ...
          patch_back_ref{plane_2,tile_2}} = [ ...
          match_points_a{patch_back_ref{plane_1,tile_1}, ...
          patch_back_ref{plane_2,tile_2}}, match_points(1,2).points_1];

        match_points_b{patch_back_ref{plane_1,tile_1}, ...
          patch_back_ref{plane_2,tile_2}} = [ ...
          match_points_b{patch_back_ref{plane_1,tile_1}, ...
          patch_back_ref{plane_2,tile_2}}, match_points(1,2).points_2];

        match_error(patch_back_ref{plane_1,tile_1}, patch_back_ref{plane_2,tile_2}) = ...
          err;
      end
    end
  end
  print_debug('---\n');
end

% make matches symmetric
for patch_id_1 = 1:n_patch_global
  for patch_id_2 = 1:patch_id_1
    if((size(match_points_a{patch_id_1, patch_id_2},2)<5 && ...
        size(match_points_b{patch_id_2, patch_id_1},2)<5) || ...
        match_error(patch_id_1, patch_id_2) > 1.0)
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

    if(~isempty(match_points_a{patch_id_1, patch_id_2}))
      print_debug('match_points_a{%d, %d}:\n', patch_id_1, patch_id_2);
      disp_debug(match_points_a{patch_id_1, patch_id_2});
      print_debug('match_points_b{%d, %d}:\n', patch_id_1, patch_id_2);
      disp_debug(match_points_b{patch_id_1, patch_id_2});
    end
    if(~isempty(match_points_a{patch_id_2, patch_id_1}))
      print_debug('match_points_a{%d, %d}:\n', patch_id_2, patch_id_1);
      disp_debug(match_points_a{patch_id_2, patch_id_1});
      print_debug('match_points_b{%d, %d}:\n', patch_id_2, patch_id_1);
      disp_debug(match_points_b{patch_id_2, patch_id_1});
    end
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
connected_component_partition = double(get_connected_components(adjacency_matrix, 4));
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

% compute affine transformations for the patches belonging to largest
% component
alpha = ones(size(match_points_1_a));

if(~isfield(sift_config, 'optimization'))
  sift_config.optimization.method = 'least_squares';
end
switch(sift_config.optimization.method)
  case 'least_squares'
    if(sift_config.is_verbose)
      fprintf('Performing Least-Squares based optimization ...');
    end
    [a,b] = getSIFTmodel(match_points_1_a, match_points_1_b, alpha, [], 1);
    if(sift_config.is_verbose)
      fprintf('done\n');
    end
  case 'levenberg_marquadt_rigid'
    if(sift_config.is_verbose)
      fprintf('Performing Levenberg-Marquadt based optimization ...');
    end
    [a,b] = ...
      get_affine_levenberg_marquadt_rigid(match_points_1_a, match_points_1_b);
    if(sift_config.is_verbose)
      fprintf('done\n');
    end
end
D = [a, b];

if(sift_config.is_verbose)
  fprintf('Collecting together the global transformations ...\n');
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
          if(sift_config.is_verbose)
            fprintf('Found a patch that did not have its own transform. Therefore, \n');
            fprintf('copying a transform from another patch in this tile.\n');
          end
          copy_patch_id_1 = original_to_new_map(patch_back_ref{plane,tile}(i(1)));
          transforms(plane,tile).transform(:,end+1) = ...
            [reshape(D{copy_patch_id_1}, [1 4]), ...
            D{copy_patch_id_1+n_member_largest_component}']';
        else
          if(sift_config.is_verbose)
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
        print_debug('ok.\n');
      end
    end
  end
end
if(sift_config.is_verbose)
  fprintf('done\n');
end

%Save files
if(sift_config.is_verbose)
  fprintf('Saving transformations ... ');
end
tform_dir = [get_region_dir(config), config.SIFT.dir];
if(~isempty(transforms))
  for plane=1:size(transforms,1)
    for tile=1:size(transforms,2)
      check_for_dir([tform_dir, image_sub_dirs{plane, tile}]);
      file_name = [tform_dir, image_prefixes{plane, tile}, ...
        '.affine_transforms.mat'];
      transform=transforms(plane, tile).transform; %#ok<NASGU>
      save2(file_name,'transform');
    end
  end
end
if(sift_config.is_verbose)
  fprintf('done\n');
end

return
end
