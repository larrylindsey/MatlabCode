function align_stack_SIFT_affine_inter_plane(config)
% align_stack_SIFT_affine_in_plane(config)
% Compute affine transformations across each plane/section based on SIFT
% correspondences. This is a two stage alignment. One affine transformation
% for each section is computed; the tranformations among tiles within each
% section are kept fixed. 
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
% v2  09192008  split SIFT affine transform into within and inter
%                 plane/section.
%

stack_config = config.stack;
sift_config = config.SIFT;

% %THINGS YOU NEED TO KNOW
% sift_config.trakEM2.getAdjacency=1;
% sift_config.makeTransforms=1;
% sift_config.makeMatches=1;
% sift_config.makeFrames=1;
% sift_config.save=1;

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
  fprintf('Computing inter plane/section transformations ...\n');
end
%%% Load the feature points, matches and within section tile transformations
%%% and transform the feature points
%%% All feature points within a section are concatenated. The SIFTtransform
%%% routine is made to think of the entire section as one "tile".
% Each z and z+1 section is made adjacent
adjacencies = ...
  diag(ones(1,size(tforms,1)-1),1) + diag(ones(1,size(tforms,1)-1),-1);

% load frames and transform them;
if(sift_config.is_verbose)
  fprintf('Loading feature points and applying in-plane transforms ...\n');
end
frames_global = {};
feature_point_id_offsets = zeros(size(tforms,1));
% plot_color = {'b', 'r', 'g', 'm', 'c'};
for plane = 1:size(tforms,1) % for each plane/section
  if(sift_config.is_verbose)
    fprintf('%d ', plane);
  end
  frames_global{plane} = [];
  current_feature_point_id_offset = 0;
  for tile = 1:size(tforms,2)
    if(isempty(tforms{plane, tile}))
      continue;
    end
    file_name=[sift_dir, image_prefixes{plane, tile}, '.sift_landmarks.mat'];
    load(file_name,'frame');

    file_name=[sift_dir, image_prefixes{plane, tile}, '.sift_transform.in_plane.mat'];
    load(file_name,'transform_p');
    transform_p = transform_p([1 3 5; 2 4 6]);
    frame(1:2, :) = transform_p * [frame(1:2, :); ones(1, size(frame,2))]; %#ok<NODEF>
%     figure(plane); plot(frame(1,:), frame(2,:), [plot_color{mod(tile,5)+1},'.']); hold on;
    frames_global{plane} = [frames_global{plane}, frame];
    feature_point_id_offsets(plane, tile) = current_feature_point_id_offset;
    current_feature_point_id_offset = current_feature_point_id_offset + size(frame,2);
  end
end
if(sift_config.is_verbose)
  fprintf('done.\n');
end

% load the matches and apply offsets on their ids
if(sift_config.is_verbose)
  fprintf('Loading matches ... \n');
end
matches_global = {};
for plane_1 = 1:size(tforms,1) % for each plane/section
  if(sift_config.is_verbose)
    fprintf('%d ', plane_1);
  end
  if(case_ids_ref(plane_1)==-1)
    continue;
  end
  for plane_2 = 1:size(tforms,1)
    matches_global{plane_1, plane_2} = [];
    if(case_ids_ref(plane_2)==-1)
      continue;
    end
    if(adjacencies(plane_1, plane_2)==0)
      continue;
    end
    file_name=[sift_dir, 'sift_matches_', num2str(case_ids_ref(plane_1)), '_', ...
      num2str(case_ids_ref(plane_2)), '.mat'];
    load(file_name,'matches');
    for tile_1 = 1:size(tforms,2)
      for tile_2 = 1:size(tforms,2)
        if(isempty(matches{tile_1, tile_2})) %#ok<USENS>
          continue;
        end
        if(size(matches{tile_1, tile_2},2)<sift_config.min_n_matches) %#ok<USENS>
          continue;
        end
        matches_global{plane_1, plane_2} = [matches_global{plane_1, plane_2}, ...
          matches{tile_1, tile_2} + ...
          [ repmat(feature_point_id_offsets(plane_1, tile_1), [1 size(matches{tile_1, tile_2},2)]); ...
            repmat(feature_point_id_offsets(plane_2, tile_2), [1 size(matches{tile_1, tile_2},2)])
          ] ];
      end
    end
  end
end
if(sift_config.is_verbose)
  fprintf('done.\n');
end

%Select best RANSAC model
if(sift_config.is_verbose)
  fprintf('Computing inter-plane transformations\n');
end
err=Inf;
if(isfield(sift_config, 'n_initial_samples'))
  options.n_initial_samples = sift_config.n_initial_samples;
end
options.p = .5;
options.w = 0.01;
options.athr = 0.01;
options.nthr = 20;
plane_offset = find(case_ids_ref~=-1, 1, 'first')-1;  % number of planes skipped from the bottom
frames_global_1 = {};
matches_global_1 = {};
images_global_1 = {};
j=0;
for i = 1:size(tforms,1)
  if(case_ids_ref(i)~=-1)
    j = j + 1;
    frames_global_1{j} = frames_global{i};
    l=0;
    for k = 1:size(tforms,1)
      if(case_ids_ref(k)~=-1)
        l = l + 1;
        matches_global_1{j,l} = matches_global{i,k};
      end
    end
    images_global_1{j} = zeros(round(max(frames_global{i}(2,:))), round(max(frames_global{i}(1,:))));
  end
end
for i=1:sift_config.n_ransac_iter
  [D1,err1] = ...
    getSIFTtransform(frames_global_1, matches_global_1, [], [], images_global_1, options);

  if(err1<err)
    D=D1;
    err=err1;
  end
  fprintf('RANSAC pass %i [error %.4g]\n',i,err1);
end
fprintf('RANSAC best pass [error %.4g]\n',err);

%Remap transforms
for plane = 1:size(tforms,1)-plane_offset % for each plane/section
  transform_g = [reshape(D{plane}, [1 4]), D{plane+length(stack_config.case_ids)}']; %#ok<NASGU>
  file_name=[sift_dir, 'section.', num2str(case_ids_ref(plane+plane_offset)), '.sift_transform.global.mat'];
  save(file_name,'transform_g');
end
if(sift_config.is_verbose)
  fprintf('\n');
end

return
end
