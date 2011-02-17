function align_stack_SIFT_affine_tile_pair_inter_plane(config)
% align_stack_SIFT_affine_tile_pair_inter_plane(config)
% Compute affine transformations between pairs of tiles of consequtive
% sections.
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
% v3  09262008  code modified from in plane alignment to inter plane
%                 tile-pair alignment. 
%

stack_config = config.stack;
sift_config = config.align.linkage_align.SIFT;

% For backward compatibility (~Dec. 15 2008)
if(isfield(sift_config, 'n_iter'))
  sift_config.RANSAC.n_iter = sift_config.n_RANSAC_iter;
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
images = {};
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
  fprintf('Computing across plane/section transformation for tiles ...\n');
end
%%%%%%%
% RANSAC parameters
%%%%%%%
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
  if(isfield(sift_config.RANSAC, 'inlier_max_relative_deviation'))
    sift_align_options.p = sift_config.RANSAC.inlier_max_relative_deviation;
  else
    sift_align_options.p = 1; % default
  end
  if(isfield(sift_config.RANSAC, 'best_inlier_carryover_frac'))
    sift_align_options.w = sift_config.RANSAC.best_inlier_carryover_frac;
  else
    sift_align_options.w = 0.15; % default
  end
  if(isfield(sift_config.RANSAC, 'max_discrepancy_for_consistency'))
    sift_align_options.nthr = sift_config.RANSAC.max_discrepancy_for_consistency;
  else
    sift_align_options.nthr = 16; % default
  end
  if(isfield(sift_config.RANSAC, 'min_n_iter'))
    sift_align_options.cthr = sift_config.RANSAC.min_n_iter;
  else
    sift_align_options.cthr = 10; % default
  end
  if(isfield(sift_config.RANSAC, 'rank_perturb_noise'))
    sift_align_options.rank_perturb_noise = sift_config.RANSAC.rank_perturb_noise;
  else
    sift_align_options.rank_perturb_noise = 2; % default
  end
end

% Generate alignment transforms between pairs of tiles from consequtive
% plane/section
for plane_1 = 1:size(tforms,1)-1 % for each plane/section
  if(case_ids_ref(plane_1)==-1)
    continue;
  end

  plane_2 = plane_1 + 1;
  if(case_ids_ref(plane_2)==-1)
    continue;
  end
  if(sift_config.is_verbose)
    fprintf('plane: %d\n', case_ids_ref(plane_1));
  end
  
  planes = [plane_1, plane_2];
  
  % load frames and matches
  frames_p = {};
  for p = 1:length(planes)
    for tile = 1:size(tforms,2)
      if(isempty(tforms{planes(p), tile}))
        continue;
      end
      file_name=[sift_dir, image_prefixes{planes(p), tile}, '.sift_landmarks.mat'];
      load2(file_name,'frame');
      frames_p{p, tile} = frame;
    end
  end

  for tile_1 = 1:size(tforms,2)
    if(isempty(tforms{plane_1, tile_1}))
      continue;
    end
    for tile_2 = 1:size(tforms,2)
      if(isempty(tforms{plane_2, tile_2}))
        continue;
      end
      if(sift_config.is_verbose)
        fprintf('--- Tile pair ---\n');
        fprintf('%d,%d\n', tile_1, tile_2);
        fprintf('%s\n%s\n', image_prefixes{plane_1, tile_1}, ...
          image_prefixes{plane_2, tile_2});
      end
      % get adjacencies between tiles
      tforms_p = {};
      tforms_p{1,1} = tforms{plane_1, tile_1};
      tforms_p{2,1} = tforms{plane_2, tile_2};
      sizes_p = {};
      sizes_p{1,1} = sizes{plane_1, tile_1};
      sizes_p{2,1} = sizes{plane_2, tile_2};
      adjacencies_p = getSIFTalpha(tforms_p, sizes_p);
      if(adjacencies_p(1,1,2,1)==0) % the tiles do not overlap
        continue;
      end
      
      frames1 = {frames_p{1, tile_1}, frames_p{2, tile_2}};
      
      matches1{1,1} = [];
      
      file_name_suffix = get_file_name_from_tuple(sift_dir, ...
        image_prefixes{plane_1, tile_1}, image_prefixes{plane_2, tile_2}, 'sm.');
      file_name = [file_name_suffix, '.mat'];
      try
        load2(file_name,'matches');
      catch
        if(sift_config.is_verbose)
          fprintf('Match file not found\n');
        end
        matches = [];
      end
      display(matches);
      matches1{1,2} = matches; %#ok<USENS>
      
      file_name_suffix = get_file_name_from_tuple(sift_dir, ...
        image_prefixes{plane_2, tile_2}, image_prefixes{plane_1, tile_1}, 'sm.');
      file_name = [file_name_suffix, '.mat'];
      try
        load2(file_name,'matches');
      catch
        if(sift_config.is_verbose)
          fprintf('Match file not found\n');
        end
        matches = [];
      end
      display(matches);
      matches1{2,1} = matches;
      
      matches1{2,2} = [];
      
      if(size(matches1{1,2},2) < sift_align_options.Nthr)
        if(sift_config.is_verbose)
          fprintf('Not enough matches\n');
          fprintf('Actual number of matches = %d, minimum number required = %d\n', ...
            size(matches1{1,2},2), sift_align_options.Nthr);
        end
        continue;
      end
      
      if(~isempty(images))
        images1 = {images{plane_1, tile_1}, images{plane_2, tile_2}};
      else
        images1 = {};
      end
      
      %Select best RANSAC model
      err=Inf;
      for i=1:sift_config.RANSAC.n_iter
        [D1, err1, f_m_p1] = ...
          getSIFTtransform(frames1,matches1,[],[],images1,sift_align_options);

        if(err1<err)
          D=D1;
          err=err1;
          final_match_points = f_m_p1;
        end
        fprintf('RANSAC pass %i [error %.4g]\n',i,err1);
      end
      fprintf('RANSAC best pass [error %.4g]\n',err);
      if(err>10 || err<0.001) % error is too high for linkage or suspiciously low
        fprintf('Suspicious error value; not saving this transformation\n');
        continue;
      end
      % abnormal change in area due to transformation - likely to be
      % erroneous
      area_factor_1 = abs(det(D{1}));
      area_factor_2 = abs(det(D{2}));
      if(min(area_factor_1,area_factor_2)<0.8)
        fprintf('Abnormal change in area; not saving this transformation\n');
        continue;
      end
      if(max(area_factor_1,area_factor_2)>1.2)
        fprintf('Abnormal change in area; not saving this transformation\n');
        continue;
      end
      
      transforms_tp = [];
      transforms_tp(1,:) = [reshape(D{1}, [1 4]), D{3}'];
      transforms_tp(2,:) = [reshape(D{2}, [1 4]), D{4}']; %#ok<NASGU>
      match_points = final_match_points; %#ok<NASGU>
      file_name_suffix = get_file_name_from_tuple(sift_dir, ...
        image_prefixes{plane_1, tile_1}, image_prefixes{plane_2, tile_2}, 'st.');
      file_name = [file_name_suffix, '.mat'];
      if(sift_config.is_verbose)
        fprintf('Saving transform and match:%s\n', file_name);
      end
      save2(file_name,'transforms_tp', 'err', 'match_points');
    end
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
