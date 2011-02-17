function align_stack_deformable_mesh_tile_pair_in_plane_SAB(config)
% align_stack_deformable_mesh_tile_pair_in_plane_SAB(config)
% Compute piecewise affine transformations using deformable mesh between
% tiles within a section. May be used for global alignment in the presence
% of folds.
% Generate stand-alone batch script.
%
% Input:
%   config    config datastructure of the reconstruction
%
% Deformable mesh code in code/aligment/tile_align_fold_deformable_mesh/
% Lou Scheffer
% Janelia Farm Research Campus, HHMI.
%
% Wrapped for pipeline by Shiv Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  08252008  init. code
% v1  08312008  modifications for pipeline - Shiv N. Vitaladevuni, JFRC, HHMI
% v2  09192008  split into SIFT feature point, matching and transform.
% v3  09262008  code modified from in plane alignment to inter plane
%                 tile-pair alignment. 
% v4  12112008  code modified for deformable mesh algorithm
% v5  12312008  modified for within section alignment
%

stack_config = config.stack;
dmesh_config = config.align.in_section_tile_pair_align.deformable_mesh;

image_dir = get_stack_dir(config);
dmesh_dir = get_deformable_mesh_dir(config);

if(~isfield(dmesh_config, 'is_verbose'))
  dmesh_config.is_verbose = true;
end
if(dmesh_config.is_verbose)
  fprintf('START: align_stack_deformable_mesh_tile_pair_in_plane_SAB\n');
  fprintf('Generating batch script for deformable mesh-based\n');
  fprintf('\twithin-plane alignment with folds ..\n');
end

% Generate alignment transforms between pairs of tiles from consequtive
% plane/section
fold_dir = get_fold_dir(config);

% Copy required binaries to the job directory
binary_names = {'deformable_mesh_cluster'};
copy_binaries_to_job_dir(config, binary_names);
% Initialize job script
fout_job = fopen_job_script(config);
fprintf(fout_job, '\n# START: align_stack_deformable_mesh_tile_pair_in_plane_SAB\n');
fprintf(fout_job, 'rm -f deformable_mesh_in_plane.long.log\n');
for plane = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(plane);
  if(dmesh_config.is_verbose)
    fprintf('case_id: %d\n', case_id);
  end
  
  fprintf(fout_job, 'echo "case_id: %d"\n', case_id);
  % get the default parameters
  params = dmesh_config.params;
  %now get any special ones for this layer
  if(isfield(dmesh_config, 'exceptions') && ...
      dmesh_config.exceptions)
    exc = dmesh_config.exceptions;
    for e=1:size(exc,1)
      if (exc{e,1}(1) == case_id && exc{e,1}(2) == case_id)
        params = [params, ' ', exc{e,2}]; %#ok<AGROW>
      end
    end
  end
  
  [image_prefixes, image_sub_dirs, is_to_be_processed] = ...
    get_image_prefixes_subdirs(config, case_id);
  tforms = get_image_tforms(config, case_id);
  sizes = get_image_sizes(config, case_id);
  
  %
  % align using deformable mesh if adjacent
  %
  for tile_1 = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_1}))
      continue;
    end
    for tile_2 = tile_1+1:length(image_prefixes)
      if(isempty(image_prefixes{tile_2}))
        continue;
      end
    
      if(dmesh_config.is_verbose)
        fprintf('%d,%d ', tile_1, tile_2);
        fprintf('--- Tile pair ---\n');
        fprintf('%s\n%s\n', image_prefixes{tile_1}, ...
          image_prefixes{tile_2});
      end
      
      if(is_to_be_processed(tile_1)==0 && is_to_be_processed(tile_2)==0)
        if(dmesh_config.is_verbose)
          fprintf('both tiles have is_to_be_processed=false, skipping this pair\n');
        end
        continue;
      end

      file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
        image_prefixes{tile_1}, image_prefixes{tile_2}, 'dt.');
      file_name = get_storage_file_name([file_name_suffix, '.mat']);
      if(dmesh_config.is_verbose)
        fprintf('Deleting existing transform file.\n');
      end
      delete(file_name);

      if(are_overlapping_affine_tforms(tforms{tile_1}, tforms{tile_2}, ...
          sizes{tile_1}, sizes{tile_2})==0)
        if(dmesh_config.is_verbose)
          fprintf('Tile pair not overlapping, skipping this pair\n');
        end
        continue;
      end
      
      fprintf(fout_job, 'echo "tile pair: %d, %d"\n', tile_1, tile_2);
      fprintf(fout_job, 'echo "%s"\n', image_prefixes{tile_1});
      fprintf(fout_job, 'echo "%s"\n', image_prefixes{tile_2});
      %repeat to add to the log file
      fprintf(fout_job, 'echo "tile pair: %d, %d" >> deformable_mesh_in_plane.long.log\n', tile_1, tile_2);
      fprintf(fout_job, 'echo "%s" >> deformable_mesh_in_plane.long.log\n', image_prefixes{tile_1});
      fprintf(fout_job, 'echo "%s" >> deformable_mesh_in_plane.long.log\n', image_prefixes{tile_2});

      if(isfield(dmesh_config, 'use_precomputed_overlap') && ...
          dmesh_config.use_precomputed_overlap)
        error('Inclusion of overlap masks in batch script not yet implemented.');
      end

      save_file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
        image_prefixes{tile_1}, image_prefixes{tile_2}, 'dt.');
      check_for_dir([dmesh_dir, image_sub_dirs{tile_1}]);
      check_for_dir([dmesh_dir, image_sub_dirs{tile_2}]);
      
      fprintf(fout_job, ...
        './deformable_mesh_cluster "%s" "%s" "%s" "%s" "%s" "%s" "%s" %s >> deformable_mesh_in_plane.long.log\n', ...
        [image_dir, image_prefixes{tile_2}, stack_config.image_suffix], ...
        [image_dir, image_prefixes{tile_1}, stack_config.image_suffix], ...
        [fold_dir, image_prefixes{tile_2}, '.fold_mask.tif'], ...
        [fold_dir, image_prefixes{tile_1}, '.fold_mask.tif'], ...
        get_storage_file_name([save_file_name_suffix, '.tforms.txt']), ...
        get_storage_file_name([save_file_name_suffix, '.map.tif']), ...
        config.job.log.deformable_mesh_in_plane, ...
        params);
    end
  end
  
  if(dmesh_config.is_verbose)
    fprintf('\n');
  end
end
fprintf(fout_job, '\n# END: align_stack_deformable_mesh_tile_pair_in_plane_SAB\n');
fclose(fout_job);

if(dmesh_config.is_verbose)
  fprintf('STOP: align_stack_deformable_mesh_tile_pair_in_plane_SAB\n');
end
return
end
