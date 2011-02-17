function correspond_tile_pair_within_section_deformable_mesh(config)
% correspond_tile_pair_within_section_deformable_mesh(config)
% Find correspondences between two images using deformable mesh. Geared
% towards images from the section.
%
% Input:
%   config    config datastructure of the reconstruction
%
% The deformable mesh algorithm was written by Louis Scheffer, Visiting
% Scientist, Janelia Farm Research Campus, HHMI.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  08252008  init. code
% v1  08312008  modifications for pipeline
% v2  09192008  split into SIFT feature point, matching and transform.
%

stack_config = config.stack;
dmesh_config = config.align.precompute.deformable_mesh_in_plane_overlap;

dmesh_dir = get_deformable_mesh_dir(config);

if(~isfield(dmesh_config, 'is_verbose'))
  dmesh_config.is_verbose = true;
end
if(~isfield(dmesh_config, 'is_verbose_figures'))
  dmesh_config.is_verbose_figures = false;
end
if(dmesh_config.is_verbose)
  fprintf('START: correspond_tile_pair_within_section_deformable_mesh\n');
end

for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  if(dmesh_config.is_verbose)
    fprintf('case_id: %d\n', case_id);
  end
  
  [images, image_prefixes, image_sub_dirs, is_to_be_processed] = ...
    get_image_from_stack(config, case_id); %#ok<ASGLU>
  
  if(isfield(dmesh_config, 'filter_version') && ~isempty(dmesh_config.filter_version))
    for j = 1:length(images)
      images{j} = filter_image2(images{j}, dmesh_config.filter_version);
    end
  end
  
  %Generate correspondences from deformable mesh based alignment
  match_points = [];
  for tile_id_1 = 1:size(images,2)
    if(isempty(images{tile_id_1}))
      continue;
    end
    for tile_id_2 = tile_id_1+1:size(images,2)
      if(isempty(images{tile_id_2}))
        continue;
      end
      
      if(dmesh_config.is_verbose)
        fprintf('%d,%d ', tile_id_1, tile_id_2);
        fprintf('--- Tile pair ---\n');
        fprintf('%s\n%s\n', image_prefixes{tile_id_1}, ...
          image_prefixes{tile_id_2});
      end
      
      if(is_to_be_processed(tile_id_1)==0 && ...
          is_to_be_processed(tile_id_2)==0)
        if(dmesh_config.is_verbose)
          fprintf('both tiles have is_to_be_processed=0, skipping this pair\n');
        end
        continue;
      end
      
      [p1, p2] = tile_pair_within_section_overlap(images{tile_id_1}, ...
        images{tile_id_2});
      if(dmesh_config.is_verbose_figures)
        display_correspondence_points(images{tile_id_1}, images{tile_id_2}, ...
          p1, p2, -1);
        title(sprintf('plane %d, tile1 %d, tile2 %d', case_id, tile_id_1, tile_id_2));
      end
      
      match_points.points_1 = p1;
      match_points.points_2 = p2;
      file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
        image_prefixes{tile_id_1}, image_prefixes{tile_id_2}, 'dmcp.');
      file_name = [file_name_suffix, '.mat'];
      save2(file_name, 'match_points');
      
      match_points.points_1 = p2;
      match_points.points_2 = p1;
      file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
        image_prefixes{tile_id_2}, image_prefixes{tile_id_1}, 'dmcp.');
      file_name = [file_name_suffix, '.mat'];
      save2(file_name, 'match_points');
    end
  end
  
  if(dmesh_config.is_verbose)
    fprintf('done.\n');
  end
end

if(dmesh_config.is_verbose)
  fprintf('STOP: correspond_tile_pair_within_section_deformable_mesh\n');
end
return
end
