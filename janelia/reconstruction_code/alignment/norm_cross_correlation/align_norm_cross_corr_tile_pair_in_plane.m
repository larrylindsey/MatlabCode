function align_norm_cross_corr_tile_pair_in_plane(config)
% align_norm_cross_corr_tile_pair_in_plane(config)
% Compute piecewise translation transformations using normalized cross
% correlation between tiles within a ection. May be used for approximation
% to more detailed alignment or as a plan B.
%
% Input:
%   config    config datastructure of the reconstruction
%
% Shiv Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  01042009  init. code
%

stack_config = config.stack;
norm_x_corr_config = config.align.precompute.norm_cross_corr_in_plane;

image_dir = get_stack_dir(config);
norm_x_corr_dir = get_norm_x_corr_dir(config);

if( ~isfield(stack_config, 'image_structure') || isempty(stack_config.image_structure) )
  error('xml file from TrakEM not specified. Could not load tile images and initial transformations.');
end;

if(~isfield(norm_x_corr_config, 'is_verbose'))
  norm_x_corr_config.is_verbose = true;
end
if(~isfield(norm_x_corr_config, 'is_verbose_figures'))
  norm_x_corr_config.is_verbose_figures = false;
end
if(norm_x_corr_config.is_verbose)
  fprintf('Normalized cross correlation-based within-plane alignment ..\n');
end

s = xmlread([image_dir, stack_config.image_structure]);
xmlwrite('~tmp.xml', s);
xmlstr = fileread('~tmp.xml');
s = xml_parseany(xmlstr);
s_t2_layer = s.t2_layer_set{1}.t2_layer;

if(norm_x_corr_config.is_verbose)
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
      dir_length = [ zeros(1, 3-length(dir_length)), dir_length]; %#ok<AGROW>

      % File path in .xml file is
      % <stack_directory>/<grid_id>/<section_id>/<image_file_name>
      image_prefixes{i,j} = s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path(dir_length(end-2)+1:end-4); %#ok<AGROW>
      image_sub_dirs{i,j} = s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path(dir_length(end-2)+1:dir_length(end)); %#ok<AGROW>

      %Extract transforms
      transform_str = s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.transform;
      if(isempty(transform_str))
        tforms{i,j}=[]; %#ok<AGROW>
      else
        matrix_str_offset = strfind(s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.transform, '(') + 1;
        matrix_str_end = strfind(s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.transform, ')') - 1;
        xform_matrix = str2num(transform_str(matrix_str_offset:matrix_str_end)); %#ok<ST2NM>
        tforms{i,j} = reshape(xform_matrix, [2 3]); %#ok<AGROW>

        width_str = s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.width;
        height_str = s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.height;
        sizes{i,j}=[str2double(width_str),str2double(height_str)]; %#ok<AGROW>
      end
    end
  end
end
if(norm_x_corr_config.is_verbose)
  fprintf('done.\n');
end

%Generate adjacency matrix from trakEM2 data
if(~isempty(tforms))
  if(norm_x_corr_config.is_verbose)
    fprintf('Computing tile adjacency matrix ..\n');
  end
  adjacencies=getSIFTalpha(tforms,sizes);
  if(norm_x_corr_config.is_verbose)
    fprintf('done.\n');
  end
else
  adjacencies=[];
end

% Generate alignment transforms between pairs of tiles from consequtive
% plane/section
fold_dir = get_fold_dir(config);
for plane = 1:size(tforms,1) % for each plane/section
  if(norm_x_corr_config.is_verbose)
    fprintf('plane: %d\n', plane);
  end
  if(case_ids_ref(plane)==-1)
    continue;
  end

  %
  % load images
  %
  images = {};
  fold_masks = {};
  for j = 1:size(tforms,2)
    if(isempty(tforms{plane, j}))
      continue;
    end
    file_name=[image_dir, image_prefixes{plane, j}, stack_config.image_suffix];
    images{j} = im2double(imread(file_name)); %#ok<AGROW>
    if(size(images{j},3)>1)
      images{j} = rgb2gray(images{j}); %#ok<AGROW>
    end
    if(max(images{j}(:))>1)
      images{j} = images{j}/max(images{j}(:)); %#ok<AGROW>
    end
    if(isfield(norm_x_corr_config, 'filter_version') && ...
        ~isempty(norm_x_corr_config.filter_version))
      images{j} = filter_image(images{j}, norm_x_corr_config.filter_version); %#ok<AGROW>
    end
    if(isfield(norm_x_corr_config, 'is_fold_considered') && ...
        norm_x_corr_config.is_fold_considered)
      file_name_prefix = [fold_dir, image_prefixes{plane, j}, '.fold_mask'];
      fold_masks{j} = load_fold_mask(file_name_prefix, config); %#ok<AGROW>
    end
  end
  
  %
  % align using normalized cross correlation if adjacent
  %
  for tile_1 = 1:length(images)
    if(isempty(images{tile_1}))
      continue;
    end
    for tile_2 = tile_1+1:length(images)
      if(isempty(images{tile_2}))
        continue;
      end
    
      if(adjacencies(plane, tile_1, plane, tile_2)==0)
        continue;
      end
      
      if(norm_x_corr_config.is_verbose)
        fprintf('%d,%d ', tile_1, tile_2);
        fprintf('--- Tile pair ---\n');
        fprintf('%s\n%s\n', image_prefixes{plane, tile_1}, ...
          image_prefixes{plane, tile_2});
      end
      if(isfield(norm_x_corr_config, 'is_fold_considered') && ...
          norm_x_corr_config.is_fold_considered)
        transform = align_translation_tile_pair_norm_cross_correlation(...
          images{tile_1}, images{tile_2}, norm_x_corr_config.scale, ...
          fold_masks{tile_1}, fold_masks{tile_2});
      else
        transform = align_translation_tile_pair_norm_cross_correlation(...
          images{tile_1}, images{tile_2}, norm_x_corr_config.scale);
      end
      if(transform.corr_value>0)
        transforms_tp = transform.transform; %#ok<NASGU>
        t = transform.transform;
        if(norm_x_corr_config.is_verbose_figures)
          display_image = zeros([size(images{tile_1})+2*size(images{tile_2}),3]);
          display_image(size(images{tile_2},1)+1:...
            size(images{tile_2},1)+size(images{tile_1},1), ...
            size(images{tile_2},2)+1:...
            size(images{tile_2},2)+size(images{tile_1},2),1) = 1-images{tile_1};
          display_image(size(images{tile_2},1)+t(6):...
            (size(images{tile_2},1)+t(6)+size(images{tile_2},1)-1), ...
            size(images{tile_2},1)+t(5):...
            (size(images{tile_2},1)+t(5)+size(images{tile_2},2)-1),2) = ...
            1-images{tile_2};
          figure; imshow(display_image);
          title(sprintf('plane %d, tile1 %d, tile2 %d\n', plane, tile_1, tile_2));
        end
        t(5:6) = -t(5:6);
        transforms_tp_rev = t;
      else
        transforms_tp = []; %#ok<NASGU>
        transforms_tp_rev = [];
      end
      file_name_suffix = get_file_name_from_tuple(norm_x_corr_dir, ...
        image_prefixes{plane, tile_1}, image_prefixes{plane, tile_2}, 'nc.');
      file_name = [file_name_suffix, '.mat'];
      save2(file_name, 'transforms_tp');
      
      transforms_tp = transforms_tp_rev; %#ok<NASGU>
      file_name_suffix = get_file_name_from_tuple(norm_x_corr_dir, ...
        image_prefixes{plane, tile_2}, image_prefixes{plane, tile_1}, 'nc.');
      file_name = [file_name_suffix, '.mat'];
      save2(file_name, 'transforms_tp');
    end
  end
  
  if(norm_x_corr_config.is_verbose)
    fprintf('\n');
  end
end

return
end
