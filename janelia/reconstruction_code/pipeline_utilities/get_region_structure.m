function region_structure = get_region_structure(config)
% region_structure = get_region_structure(config)
%
% Get the structure of the region.
% config.region.region_structure: contains all the images in the region.
%   .planes(1:n_plane): planes(p) has info. regarding plane pth in region.
%     .z: z depth (double)
%     .thickness: section thickness (double)
%     .tilt_planes(:): tilt_planes(tp) has info for tp tilt-plane in planes(p).
%       .tilt_angle:  tilt angle for tilt-plane tp.
%       .tiles(:): tiles(t) has info for t tile in tilt-plane(tp)
%         .name:  image prefix. E.g, if full path is
%                   "<stack_dir>/a/b/c.tif" then name is "a/b/c"
%         .sub_dir: a/b/
%         .width
%         .height
%         .transform
%         .style
%         .links
%         .type
%         .o_width
%         .o_height
%         .file_path
%         .is_to_be_processed
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  10082009  init code
%

global config_global

if(config.is_verbose)
  fprintf('START: get_region_structure\n');
end

region_structure = [];
region_structure.planes = [];
stack_config = config.stack;
if(isfield(stack_config, 'image_structure') && ...
    ~isempty(stack_config.image_structure))
  file_name = [get_stack_dir(config), stack_config.image_structure];
  if(config.is_verbose)
    fprintf('Reading in %s\n', file_name);
  end
  
  % get the differential stack structure if specified
  s_diff = [];
  if(isfield(stack_config, 'image_structure_differential') && ...
      ~isempty(stack_config.image_structure_differential))
    file_name_diff = [get_stack_dir(config), ...
      stack_config.image_structure_differential];
    if(~isfield(stack_config, 'image_structure_is_kludged') || ...
        ~stack_config.image_structure_is_kludged)
      s = xmlread(file_name_diff);
      xmlwrite([config_global.temp_dir, '~tmp.xml'], s);
      xmlstr = fileread([config_global.temp_dir, '~tmp.xml']);
    else
      xmlstr = fileread(file_name_diff);
    end
    s = xml_parseany(xmlstr);
    s_t2_layer = s.t2_layer_set{1}.t2_layer;
    
    s_diff.plane_z = [];
    s_diff.tile_names = {};
    for i = 1:length(s_t2_layer)
      z = str2double(s_t2_layer{i}.ATTRIBUTE.z);
      if(ismember(z, stack_config.case_ids))
        s_diff.plane_z(end+1) = z;
        
        plane_id = length(s_diff.plane_z);
        s_diff.tile_names{plane_id} = {};
        for j = 1:length(s_t2_layer{i}.t2_patch)
          % remove the directory and extension from the file_path and store
          % in image_prefix
          dir_length = ...
            strfind(s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path, '/');
          dir_length = ...
            [zeros(1, max(0, 3-length(dir_length))), dir_length]; %#ok<AGROW>
          
          % File path in .xml file is
          % <stack_directory>/<grid_id>/<section_id>/<image_file_name>
          s_diff.tile_names{plane_id}{end+1} = ...
            s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path(dir_length(end-2)+1:end-4);
        end
      end
    end
  end
  
  if(~isfield(stack_config, 'image_structure_is_kludged') || ...
      ~stack_config.image_structure_is_kludged)
    s = xmlread(file_name);
    xmlwrite([config_global.temp_dir, '~tmp.xml'], s);
    xmlstr = fileread([config_global.temp_dir, '~tmp.xml']);
  else
    xmlstr = fileread(file_name);
  end
  s = xml_parseany(xmlstr);
  s_t2_layer = s.t2_layer_set{1}.t2_layer;
  
  for i = 1:length(s_t2_layer)
    z = str2double(s_t2_layer{i}.ATTRIBUTE.z);
    if(ismember(z, stack_config.case_ids))
      plane = [];
      plane.z = z;
      plane.thickness = ...
        str2double(s_t2_layer{i}.ATTRIBUTE.thickness);
      plane.tilt_planes = [];
      
      % for now just the 0 tilt-plane
      tilt_plane = [];
      tilt_plane.tilt_angle = 0.0;
      tilt_plane.tiles = [];
      
      for j = 1:length(s_t2_layer{i}.t2_patch)
        tile = [];
        % remove the directory and extension from the file_path and store
        % in image_prefix
        dir_length = ...
          strfind(s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path, '/');
        dir_length = ...
          [zeros(1, max(0, 3-length(dir_length))), dir_length]; %#ok<AGROW>
        
        % File path in .xml file is
        % <stack_directory>/<grid_id>/<section_id>/<image_file_name>
        tile.name = ...
          s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path(dir_length(end-2)+1:end-4);
        tile.sub_dir = ...
          s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.file_path(dir_length(end-2)+1:dir_length(end));

        tile.width = str2double(s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.width);
        tile.height = str2double(s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.height);
        
        %Extract transforms
        transform_str = s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.transform;
        if(isempty(transform_str))
          tile.transform = [];
        else
          matrix_str_offset = strfind(s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.transform, '(') + 1;
          matrix_str_end = strfind(s_t2_layer{i}.t2_patch{j}.ATTRIBUTE.transform, ')') - 1;
          xform_matrix = str2num(transform_str(matrix_str_offset:matrix_str_end)); %#ok<ST2NM>
          tile.transform = reshape(xform_matrix, [2 3]);
        end
        
        tile.is_to_be_processed = true;
        if(~isempty(s_diff))
          plane_id = find(s_diff.plane_z==z);
          if(~isempty(plane_id))
            is_in_differential = max(strcmp(tile.name, ...
              s_diff.tile_names{plane_id}))>0;
            if(is_in_differential)
              tile.is_to_be_processed = false;
            end
          end
        end
        
        if(isempty(tilt_plane.tiles))
          tilt_plane.tiles = tile;
        else
          tilt_plane.tiles(end+1) = tile;
        end
      end
      
      if(isempty(plane.tilt_planes))
        plane.tilt_planes = tilt_plane;
      else
        plane.tilt_planes(end+1) = tilt_plane;
      end
      
      if(isempty(region_structure.planes))
        region_structure.planes = plane;
      else
        region_structure.planes(end+1) = plane;
      end
    end
  end
else
  for i = 1:length(stack_config.case_ids)
    plane = [];
    plane.z = stack_config.case_ids(i);
    plane.thickness = 1.0; %default thickness
    plane.tilt_planes = [];
    
    % for now just the 0 tilt-plane
    tilt_plane = [];
    tilt_plane.tilt_angle = 0.0;
    tilt_plane.tiles = [];
    
    tile.name = sprintf(stack_config.image_prefix, plane.z);
    tile.sub_dir = '';
    tile.width = -1; % for now invalidate
    tile.height = -1; % for now invalidate
    tile.transform = [1 0 0; 0 1 0];  % identity affine transform
    tile.is_to_be_processed = true;
    
    if(isempty(tilt_plane.tiles))
      tilt_plane.tiles = tile;
    else
      tilt_plane.tiles(end+1) = tile;
    end
    
    if(isempty(plane.tilt_planes))
      plane.tilt_planes = tilt_plane;
    else
      plane.tilt_planes(end+1) = tilt_plane;
    end
    
    if(isempty(region_structure.planes))
      region_structure.planes = plane;
    else
      region_structure.planes(end+1) = plane;
    end
  end
end

if(config.is_verbose)
  fprintf('STOP: get_region_structure\n');
end

return;
end
