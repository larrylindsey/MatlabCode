function align_stack_deformable_mesh_tile_pair_inter_plane_SAB(config)
% align_stack_deformable_mesh_tile_pair_inter_plane_SAB(config)
% Compute piecewise affine transformations using deformable mesh between
% pairs of tiles of consequtive sections.
% Generate stand alone batch script.
%
% Input:
%   config    config datastructure of the reconstruction
%
% Deformable mesh code in code/aligment/tile_align_fold_deformable_mesh/
% Lou Scheffer
% Visiting Scientist, Janelia Farm Research Campus, HHMI.
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
%

stack_config = config.stack;
dmesh_config = config.align.linkage_align.deformable_mesh;

image_dir = get_stack_dir(config);
dmesh_dir = get_deformable_mesh_dir(config);

if(~isfield(dmesh_config, 'is_verbose'))
  dmesh_config.is_verbose = true;
end
if(dmesh_config.is_verbose)
  fprintf('START: align_stack_deformable_mesh_tile_pair_inter_plane_SAB\n');
end

% Generate alignment transforms between pairs of tiles from consequtive
% plane/section
fold_dir = get_fold_dir(config);

% Copy required binaries to the job directory
binary_names = {'deformable_mesh_cluster'};
copy_binaries_to_job_dir(config, binary_names);
% Initialize job script
fout_job = fopen_job_script(config);
fprintf(fout_job, '\n# START: align_stack_deformable_mesh_tile_pair_inter_plane_SAB\n');
fprintf(fout_job, 'rm -f deformable_mesh_inter_plane.long.log\n');
for plane_1 = 1:length(stack_config.case_ids)-1 % for each plane/section
  case_id_1 = stack_config.case_ids(plane_1);
  if(dmesh_config.is_verbose)
    fprintf('case_id_1: %d\n', case_id_1);
  end
  
  fprintf(fout_job, 'echo "plane_1: %d"\n', case_id_1);

  case_id_2 = case_id_1 + 1;
  if(~ismember(case_id_2, stack_config.case_ids))
    continue;
  end
  fprintf(fout_job, 'echo "plane_2: %d"\n', case_id_2);
  % get the default parameters
  % params1 will be used for the first of the two calls, params2 for the
  % second.  First they both get the general parameters
  params1 = dmesh_config.params;
  params2 = params1;
  %now get any special ones for this layer pair
  if(isfield(dmesh_config, 'exceptions') && ...
      ~isempty(dmesh_config.exceptions))
    exc = dmesh_config.exceptions;
    for e=1:size(exc,1)
      if (exc{e,1}(1) == case_id_1 && exc{e,1}(2) == case_id_2)
        params1 = [params1, ' ', exc{e,2}]; %#ok<AGROW>
      end
      if (exc{e,1}(1) == case_id_2 && exc{e,1}(2) == case_id_1)
        params2 = [params2, ' ', exc{e,2}]; %#ok<AGROW>
      end
    end
    fprintf(fout_job, 'echo params1="%s" params2="%s"\n', params1, params2);
  end
  
  [image_prefixes_1, image_sub_dirs_1, is_to_be_processed_1] = ...
    get_image_prefixes_subdirs(config, case_id_1);
  tforms_1 = get_image_tforms(config, case_id_1);
  sizes_1 = get_image_sizes(config, case_id_1);

  [image_prefixes_2, image_sub_dirs_2, is_to_be_processed_2] = ...
    get_image_prefixes_subdirs(config, case_id_2);
  tforms_2 = get_image_tforms(config, case_id_2);
  sizes_2 = get_image_sizes(config, case_id_2);
  %
  % align using deformable mesh if adjacent
  %
  for tile_1 = 1:length(image_prefixes_1)
    if(isempty(image_prefixes_1{tile_1}))
      continue;
    end
    for tile_2 = 1:length(image_prefixes_2)
      if(isempty(image_prefixes_2{tile_2}))
        continue;
      end
    
      if(dmesh_config.is_verbose)
        fprintf('--- Tile pair ---\n');
        fprintf('%d,%d\n', tile_1, tile_2);
        fprintf('%s\n%s\n', image_prefixes_1{tile_1}, ...
          image_prefixes_2{tile_2});
      end
      
      if(is_to_be_processed_1(tile_1)==0 && is_to_be_processed_2(tile_2)==0)
        if(dmesh_config.is_verbose)
          fprintf('both tiles have is_to_be_processed=false, skipping this pair\n');
        end
        continue;
      end
      
      if(dmesh_config.is_verbose)
        fprintf('Deleting pre-existing transform file.\n');
      end
      file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
        image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'dt.');
      file_name = get_storage_file_name([file_name_suffix, '.mat']);
      delete(file_name);
      file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
        image_prefixes_2{tile_2}, image_prefixes_1{tile_1}, 'dt.');
      file_name = get_storage_file_name([file_name_suffix, '.mat']);
      delete(file_name);
      
      if(are_overlapping_affine_tforms(tforms_1{tile_1}, tforms_2{tile_2}, ...
          sizes_1{tile_1}, sizes_2{tile_2})==0)
        if(dmesh_config.is_verbose)
          fprintf('tiles are not overlapping, skipping this pair\n');
        end
        continue;
      end
      
      fprintf(fout_job, 'echo "tile pair: %d, %d"\n', tile_1, tile_2);
      fprintf(fout_job, 'echo "%s"\n', image_prefixes_1{tile_1});
      fprintf(fout_job, 'echo "%s"\n', image_prefixes_2{tile_2});
      %repeat, adding to the log file
      fprintf(fout_job, 'echo "tile pair: %d, %d" >> deformable_mesh_inter_plane.long.log\n', tile_1, tile_2);
      fprintf(fout_job, 'echo "%s" >> deformable_mesh_inter_plane.long.log\n', image_prefixes_1{tile_1});
      fprintf(fout_job, 'echo "%s" >> deformable_mesh_inter_plane.long.log\n', image_prefixes_2{tile_2});

      if(isfield(dmesh_config, 'use_precomputed_overlap') && ...
          dmesh_config.use_precomputed_overlap)
        error('Inclusion of overlap masks in batch script not yet implemented.');
      end
      
      % flip center angle for deformable mesh program if the sections are
      % rotated by 180 degrees
      if(tforms_1{tile_1}(1)*tforms_2{tile_2}(1)>=0)
        center_angle = 0;
      else
        center_angle = 180;
      end
      
      save_file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
        image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'dt.');
      check_for_dir([dmesh_dir, image_sub_dirs_1{tile_1}]);
      check_for_dir([dmesh_dir, image_sub_dirs_2{tile_2}]);

      fprintf(fout_job, ...
        './deformable_mesh_cluster "%s" "%s" "%s" "%s" "%s" "%s" "%s" -CTR=%d %s >> deformable_mesh_inter_plane.long.log\n', ...
        [image_dir, image_prefixes_2{tile_2}, stack_config.image_suffix], ...
        [image_dir, image_prefixes_1{tile_1}, stack_config.image_suffix], ...
        [fold_dir, image_prefixes_2{tile_2}, '.fold_mask.tif'], ...
        [fold_dir, image_prefixes_1{tile_1}, '.fold_mask.tif'], ...
        get_storage_file_name([save_file_name_suffix, '.tforms.txt']), ...
        get_storage_file_name([save_file_name_suffix, '.map.tif']), ...
        config.job.log.deformable_mesh_inter_plane, ...
        center_angle, ...
        params1);

      save_file_name_suffix = get_file_name_from_tuple(dmesh_dir, ...
        image_prefixes_2{tile_2}, image_prefixes_1{tile_1}, 'dt.');

      fprintf(fout_job, ...
        './deformable_mesh_cluster "%s" "%s" "%s" "%s" "%s" "%s" "%s" -CTR=%d %s >> deformable_mesh_inter_plane.long.log\n', ...
        [image_dir, image_prefixes_1{tile_1}, stack_config.image_suffix], ...
        [image_dir, image_prefixes_2{tile_2}, stack_config.image_suffix], ...
        [fold_dir, image_prefixes_1{tile_1}, '.fold_mask.tif'], ...
        [fold_dir, image_prefixes_2{tile_2}, '.fold_mask.tif'], ...
        get_storage_file_name([save_file_name_suffix, '.tforms.txt']), ...
        get_storage_file_name([save_file_name_suffix, '.map.tif']), ...
        config.job.log.deformable_mesh_inter_plane, ...
        center_angle, ...
        params2);
    end
  end
  
  if(dmesh_config.is_verbose)
    fprintf('\n');
  end
end
fprintf(fout_job, '\n# END: align_stack_deformable_mesh_tile_pair_inter_plane_SAB\n');
fclose(fout_job);

if(dmesh_config.is_verbose)
  fprintf('STOP: align_stack_deformable_mesh_tile_pair_inter_plane_SAB\n');
end

return
end
