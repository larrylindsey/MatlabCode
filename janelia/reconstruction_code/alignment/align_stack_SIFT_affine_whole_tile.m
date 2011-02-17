function transforms = align_stack_SIFT_affine_whole_tile(config)
% [frames, descriptors, matches, adjacencies, transformations] = ...
%   align_stack_SIFT_affine_whole_tile(config)
% Constructs SIFT landmarks and matches for image from the stack either
% from the reconstruction's "config" datastructure or from Leginon-format
% .xml file. Compute affine transformations based on SIFT correspondences.
%
% THIS USES getSIFT* LIBRARY
%
% Input:
%   config    config datastructure of the reconstruction
% Output:
%   frames              cell array of SIFT frames (feature point
%                         coordinates)
%   descriptors         cell array of descriptors
%   matches             cell array of matches
%   adjacencies         matrix specifying image-adjacencies
%   transformations     global alignment transformations
%
% Yuriy Mishchenko
% Janelia Farm Research Campus, HHMI.
%
% v0  08252008  init. code
% v1  08312008  modifications for pipeline
%

stack_config = config.stack;
sift_config = config.align.global_align.SIFT;

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
      dir_length = [zeros(1, 3-length(dir_length)), dir_length];

      % File path in .xml file is
      % <stack_directory>/<grid_id>/<section_id>/<image_file_name>
%       images{i,j} = im2double(imread(...
%         [image_dir, s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path(dir_length(end-2)+1:end)]));
%       if(max(images{i,j}(:))>1)
%         images{i,j} = images{i,j}/max(images{i,j}(:));
%       end
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
  mkdir2(sift_dir);
  if(~isempty(adjacencies))
    file_name=[sift_dir, sprintf('sift_adjacencies.mat')];
    save2(file_name,'adjacencies');
  end
end


% load SIFT frames
frames={};
for i=1:size(tforms,1)
  for j=1:size(tforms,2)
    if(isempty(tforms{i,j}))
      continue;
    end
    file_name=[sift_dir, image_prefixes{i,j}, '.sift_landmarks.mat'];
    load2(file_name,'frame','descriptor');
    frames{i,j} = frame;
  end
end

% load SIFT matches
matches_global=cell(size(tforms,1), size(tforms,2), size(tforms,1), size(tforms,2));
for plane_1 = 1:size(tforms,1) % for each plane/section
  if(sift_config.is_verbose)
    fprintf('%d ', plane_1);
  end
  for plane_2 = 1:size(tforms,1)
    if(sum(sum(squeeze(adjacencies(plane_1, :, plane_2, :))))==0)
      continue;
    end
    for tile_1 = 1:size(tforms,2)
      if(isempty(image_prefixes{plane_1, tile_1}))
        continue;
      end
      for tile_2 = 1:size(tforms,2)
        if(isempty(image_prefixes{plane_2, tile_2}))
          continue;
        end
        file_name_suffix = get_file_name_from_tuple(sift_dir, ...
          image_prefixes{plane_1, tile_1}, image_prefixes{plane_2, tile_2}, 'sm.');
        file_name = [file_name_suffix, '.mat'];
        try
          load2(file_name,'matches');
        catch
          matches = [];
        end
        matches_global{plane_1, tile_1, plane_2, tile_2} = matches;
      end
    end
  end
end

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
%Generate global alignment transforms
transforms={};
if(~isempty(matches_global))
  if(sift_config.is_verbose)
    fprintf('Computing SIFT based affine transformations ..\n');
  end
  [frames1,matches1,ind1,adjacencies1, images1] = ...
    getSIFTrow(frames,matches_global,adjacencies,images);

  %Select best RANSAC model
  err=Inf;
  for i=1:sift_config.RANSAC.n_iter
    [D1,err1]=getSIFTtransform(frames1,matches1,ind1,adjacencies1,images1,sift_align_options);

    if(err1<err)
      D=D1;
      err=err1;
    end
    fprintf('RANSAC pass %i [error %.4g]\n',i,err1);
  end
  fprintf('RANSAC best pass [error %.4g]\n',err);

  %Remap transforms
  for plane=1:size(tforms,1)
    for tile=1:size(tforms,2)
      cnt=ind1(plane,tile);
      if(cnt==0)
        continue;
      end
      transforms{plane,tile}=[reshape(D{cnt}, [1 4]), D{cnt+max(ind1(:))}'];
    end
  end
  if(sift_config.is_verbose)
    fprintf('done.\n');
  end
end

%Save files
tform_dir = [get_region_dir(config), config.SIFT.dir];
if(~isempty(transforms))
  for i=1:size(transforms,1)
    for j=1:size(transforms,2)
      if(isempty(transforms{i,j}))
        continue;
      end
      check_for_dir([tform_dir, image_sub_dirs{i,j}]);
      file_name = [tform_dir, image_prefixes{i,j}, '.sift_transforms.mat'];
      transform=transforms{i,j}; %#ok<NASGU>
      save2(file_name,'transform');
    end
  end
end

return
end
