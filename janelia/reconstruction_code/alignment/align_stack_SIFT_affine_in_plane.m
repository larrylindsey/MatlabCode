function align_stack_SIFT_affine_in_plane(config)
% align_stack_SIFT_affine_in_plane(config)
% Compute affine transformations within each plane/section based on SIFT
% correspondences.
%
% THIS USES getSIFT* LIBRARY
%
% Input:
%   config    config datastructure of the reconstruction
%
% Yuriy Mishchenko
% Janelia Farm Research Campus, HHMI.
%
% v0  08252008  init. code
% v1  08312008  modifications for pipeline - Shiv N. Vitaladevuni, JFRC, HHMI
% v2  09192008  split into SIFT feature point, matching and transform.
%

stack_config = config.stack;
sift_config = config.align.in_section_align.SIFT;

% For backward compatibility (~Dec. 15 2008)
if(isfield(sift_config, 'n_ransac_iter'))
  sift_config.RANSAC.n_iter = sift_config.n_ransac_iter;
end
if(isfield(sift_config, 'n_initial_samples'))
  sift_config.RANSAC.n_initial_samples = sift_config.n_initial_samples;
end
if(isfield(sift_config, 'min_n_matches'))
  sift_config.RANSAC.min_n_matches = sift_config.min_n_matches;
end

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
  fprintf('SIFT-based alignment ..\n');
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
case_ids_ref = -ones(1,length(s_t2_layer)); % record of case_id of each plane to be used when saving
for i = 1:length(s_t2_layer)
  case_id=str2double(s_t2_layer{i}.ATTRIBUTE.z);
  if(ismember(case_id,stack_config.case_ids) || isempty(stack_config.case_ids))
    case_ids_ref(i) = case_id;
    for j = 1:length(s_t2_layer{i}.t2_patch)
      % remove the directory and extension from the file_path and store
      % in image_prefix
      dir_length = strfind(s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path, '/');
      dir_length = [ zeros(1, 3-length(dir_length)), dir_length];

      % File path in .xml file is
      % <stack_directory>/<grid_id>/<section_id>/<image_file_name>
      if(sift_config.is_verbose_figures)
        images{i,j} = im2double(imread(...
          [image_dir, s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path(dir_length(end-2)+1:end)]));
        if(size(images{i,j},3)>1)
          images{i,j} = rgb2gray(images{i,j});
        end
        if(max(images{i,j}(:))>1 || max(images{i,j}(:))<0.1)
          images{i,j} = images{i,j}/(max(images{i,j}(:))+eps);
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
if(sift_config.is_verbose)
  fprintf('done.\n');
end

if(sift_config.is_verbose)
  fprintf('Computing within plane/section transformation for tiles ...\n');
end
%Generate within plane/section alignment transforms
for plane = 1:size(tforms,1) % for each plane/section
  if(case_ids_ref(plane)==-1)
    continue;
  end
  if(sift_config.is_verbose)
    fprintf('plane: %d\n', case_ids_ref(plane));
  end

  % get number of tiles in section, if 1 then simply write identity
  % transformation.
  n_tile = 0;
  for tile = 1:size(tforms,2)
    if(isempty(tforms{plane, tile}))
      continue;
    end
    n_tile = n_tile + 1;
  end
  if(n_tile==1)
    image_set_string = '.';
    for tile = 1:size(tforms,2)
      if(isempty(tforms{plane, tile}))
        continue;
      end
      image_set_string = [image_set_string, image_prefixes{plane, tile}]; %#ok<AGROW>
    end
    image_set_string = strrep(image_set_string, '/', '_');
    image_set_string = strrep(image_set_string, '\', '_');
    transform_p = [1 0 0 1 0 0]; %#ok<NASGU>
    for tile = 1:size(tforms,2)
      if(isempty(tforms{plane, tile}))
        continue;
      end
      check_for_dir([sift_dir, image_sub_dirs{plane, tile}]);
      file_name=[sift_dir, image_prefixes{plane, tile}, image_set_string, ...
        '.sift_transform.in_plane.mat'];
      save2(file_name,'transform_p');
    end
    continue;
  end
  % get adjacencies between tiles
  tforms_p = {};
  sizes_p = {};
  for tile = 1:size(tforms,2)
    tforms_p{1, tile} = tforms{plane, tile};
    sizes_p{1, tile} = sizes{plane, tile};
  end
  adjacencies_p = getSIFTalpha(tforms_p,sizes_p);
  
  % load frames and matches
  frames_p = {};
  for tile = 1:size(tforms,2)
    if(isempty(tforms{plane, tile}))
      continue;
    end
    file_name=[sift_dir, image_prefixes{plane, tile}, '.sift_landmarks.mat'];
    load2(file_name,'frame');
    frames_p{1, tile} = frame;
  end
  
  matches_p = {};
  for tile_1 = 1:size(tforms,2)
    if(isempty(tforms{plane, tile_1}))
      continue;
    end
    for tile_2 = 1:size(tforms,2)
      if(isempty(tforms{plane, tile_2}))
        continue;
      end
      file_name_suffix = get_file_name_from_tuple(sift_dir, ...
        image_prefixes{plane, tile_1}, image_prefixes{plane, tile_2}, 'sm.');
      file_name = [file_name_suffix, '.mat'];
      try
        load2(file_name,'matches');
      catch
        matches = [];
      end
      matches_p{1,tile_1, 1, tile_2} = matches;
      if(sift_config.is_verbose_figures)
        if(~isempty(matches))
          display_correspondence_points(images{plane, tile_1}, ...
            images{plane, tile_2}, frames_p{1, tile_1}(1:2,matches(1,:)), ...
            frames_p{1, tile_2}(1:2,matches(2,:)), -1);
          title(sprintf('plane: %d, tile 1: %d, tile 2: %d', case_ids_ref(plane), ...
            tile_1, tile_2));
        end
      end
    end
  end

  % Rearrange the cell arrays for SIFT transform routine
  [frames1, matches1, ind1, adjacencies1, images1] = ...
    getSIFTrow(frames_p, matches_p, adjacencies_p, {});
  
  %%%%%%%
  % RANSAC parameters
  %%%%%%%
  sift_align_options = [];
  if(isfield(sift_config, 'RANSAC'))
    if(isfield(sift_config.RANSAC, 'n_initial_samples'))
      sift_align_options.n_initial_samples = sift_config.RANSAC.n_initial_samples;
    end
    if(isfield(sift_config.RANSAC, 'min_n_matches'))
      sift_align_options.Nthr = sift_config.RANSAC.min_n_matches;
    end
    if(isfield(sift_config.RANSAC, 'consistency'))
      sift_align_options.athr = sift_config.RANSAC.consistency;
    end
    if(isfield(sift_config.RANSAC, 'best_inlier_carryover_frac'))
      sift_align_options.w = sift_config.RANSAC.best_inlier_carryover_frac;
    else
      sift_align_options.w = 0.25; % default
    end
  end
  %Select best RANSAC model
  err=Inf;
  deformable_mesh_matches = [];
  if(isfield(sift_config, 'plan_B_method') && ...
      strcmp(sift_config.plan_B_method, 'deformable_mesh')==1)
    dmesh_dir = get_deformable_mesh_dir(config);
    for tile_1 = 1:size(tforms,2)
      if(isempty(image_prefixes{plane, tile_1}))
        continue;
      end
      for tile_2 = 1:size(tforms,2)
        if(isempty(image_prefixes{plane, tile_2}))
          continue;
        end
        file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
          image_prefixes{plane, tile_1}, image_prefixes{plane, tile_2}, 'dmcp.');
        file_name = [file_name_suffix, '.mat'];
        try
          load2(file_name, 'match_points');
        catch
          match_points.points_1 = [];
          match_points.points_2 = [];
        end
        deformable_mesh_matches.match_points(tile_1, tile_2) = match_points;
      end
    end
  else
    deformable_mesh_matches.match_points = [];
  end
  for i=1:sift_config.RANSAC.n_iter
    [D1,err1]=getSIFTtransform(frames1,matches1,ind1,adjacencies1,images1, ...
      sift_align_options, deformable_mesh_matches.match_points);

    if(err1<err)
      D=D1;
      err=err1;
    end
    fprintf('RANSAC pass %i [error %.4g]\n',i,err1);
  end
  fprintf('RANSAC best pass [error %.4g]\n',err);

  % Generate a substring for the file name that accounts for the images
  % involved in the joint alignment
  image_set_string = '.';
  for tile = 1:size(tforms,2)
    if(isempty(tforms{plane, tile}))
      continue;
    end
    image_set_string = [image_set_string, image_prefixes{plane, tile}];
  end
  image_set_string = strrep(image_set_string, '/', '_');
  image_set_string = strrep(image_set_string, '\', '_');

  for tile = 1:size(tforms,2)
    if(isempty(tforms{plane, tile}))
      continue;
    end
    cnt=ind1(tile);
    if(cnt==0)
      continue;
    end
    transform_p = [reshape(D{cnt}, [1 4]), D{cnt+max(ind1(:))}']; %#ok<NASGU>
    check_for_dir([sift_dir, image_sub_dirs{plane, tile}]);
    file_name=[sift_dir, image_prefixes{plane, tile}, image_set_string, ...
      '.sift_transform.in_plane.mat'];
    save2(file_name,'transform_p');
  end
  if(sift_config.is_verbose)
    fprintf('\n');
  end
end
if(sift_config.is_verbose)
  fprintf('done.\n');
end

return
end
