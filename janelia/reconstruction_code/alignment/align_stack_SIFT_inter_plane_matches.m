function matches = align_stack_SIFT_inter_plane_matches(config)
% matches = align_stack_SIFT_inter_plane_matches(config)
% Constructs SIFT inter-plane matches for image from the stack either
% from the reconstruction's "config" datastructure or from Leginon-format
% .xml file.
%
% THIS USES getSIFT* LIBRARY
%
% Input:
%   config    config datastructure of the reconstruction
% Output:
%   matches             cell array of matches
%
% Yuriy Mishchenko
% Janelia Farm Research Campus, HHMI.
%
% v0  08252008  init. code
% v1  08312008  modifications for pipeline
% v2  09192008  split into SIFT feature point, matching and transform.
%

fprintf('START: align_stack_SIFT_inter_plane_matches\n');
stack_config = config.stack;
if(strcmp(config.align.linkage_align.method, 'SIFT')~=1)
  warning('Linkage align method is not SIFT. Doing nothing.\n'); %#ok<WNTAG>
  fprintf('STOP: align_stack_SIFT_inter_plane_matches\n');
  return;
end

sift_config = config.align.linkage_align.SIFT;

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


%Generate adjacency matrix from trakEM2 data
if(~isempty(tforms))
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

if(isfield(sift_config, 'save') && sift_config.save)
  mkdir2(sift_dir);
  if(~isempty(adjacencies))
    file_name=[sift_dir, sprintf('sift_adjacencies.mat')];
    save2(file_name,'adjacencies');
  end
end


%Generate SIFT matches
frames={};
descriptors={};
for i=1:size(tforms,1)
  for j=1:size(tforms,2)
    if(isempty(tforms{i,j}))
      continue;
    end
    file_name=[sift_dir, image_prefixes{i,j}, '.sift_landmarks.mat'];
    load2(file_name,'frame','descriptor');
    frames{i,j}=frame; descriptors{i,j}=descriptor;
  end
end
if(sift_config.is_verbose)
  fprintf('Computing SIFT matches ..\n');
end
matches_all = getSIFTmatches_inter_plane(frames, descriptors);
if(sift_config.is_verbose)
  fprintf('done.\n');
end

if(~isempty(matches_all))
  for plane_1 = 1:size(matches_all, 1)
    for plane_2 = 1:size(matches_all, 3)
      for tile_1 = 1:size(matches_all, 2)
        for tile_2 = 1:size(matches_all, 4)
          matches = matches_all{plane_1, tile_1, plane_2, tile_2};
          if(~isempty(matches))
            if(sift_config.is_verbose_figures)
              im1 = imread([get_stack_dir(config), image_prefixes{plane_1, tile_1}, ...
                '.tif']);
              im2 = imread([get_stack_dir(config), image_prefixes{plane_2, tile_2}, ...
                '.tif']);
              if(isfield(config.align.precompute.SIFT.feature_point, 'filter_version') && ...
                  ~isempty(config.align.precompute.SIFT.feature_point.filter_version))
                im1 = filter_image2(im1, config.align.precompute.SIFT.feature_point.filter_version);
                im2 = filter_image2(im2, config.align.precompute.SIFT.feature_point.filter_version);
              end
              disp_scale = 3;
              im1 = imresize(im1, 1/disp_scale);
              im2 = imresize(im2, 1/disp_scale);
              disp_im = [im1, im2];
              figure;
              imshow(disp_im);
              hold on;
              plot(frames{plane_1, tile_1}(1,matches(1,:))/disp_scale, frames{plane_1, tile_1}(2,matches(1,:))/disp_scale, 'go');
              plot(frames{plane_2, tile_2}(1,matches(2,:))/disp_scale + size(im1,2), frames{plane_2, tile_2}(2,matches(2,:))/disp_scale, 'go');
              plot([frames{plane_1, tile_1}(1,matches(1,:))/disp_scale; frames{plane_2, tile_2}(1,matches(2,:))/disp_scale + size(im1,2)], ...
                [frames{plane_1, tile_1}(2,matches(1,:))/disp_scale; frames{plane_2, tile_2}(2,matches(2,:))/disp_scale], 'r');
              hold off
            end
            file_name_suffix = get_file_name_from_tuple(sift_dir, ...
              image_prefixes{plane_1, tile_1}, image_prefixes{plane_2, tile_2}, 'sm.');
            file_name = [file_name_suffix, '.mat'];
            save2(file_name,'matches');
          end
        end
      end
    end
  end
end

fprintf('STOP: align_stack_SIFT_inter_plane_matches\n');
return
end
