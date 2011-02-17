function dump_info_for_c_rendering(config)
% dump_info_for_c_rendering(config)
% Dump data for c-based rendering implemented by Lou.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
% v1  09302008  modified to multiple tile_id per section (trakEM xml)
%

if(config.is_verbose)
  fprintf('START: generate_matlab_global_tforms_multi_tile_patches\n');
end

stack_config = config.stack;
switch(config.align.global_align.method)
  case 'SIFT'
    error('This is not meant for SIFT-based alignment');
  case 'deformable_mesh'
    tform_dir = [get_region_dir(config), config.deformable_mesh.dir];
    tform_save_dir = [get_to_be_proofread_dir(config), ...
      config.deformable_mesh.dir];
end

check_for_dir(get_global_stitching_param_dir(config));

fout_info = fopen([get_region_dir(config), 'render_info.txt'], 'wt');
for layer_id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(layer_id);
  fprintf('plane: %d\n', case_id);
  
  [image_prefixes, image_sub_dirs]  = ...
    get_image_prefixes_subdirs(config, case_id);
  
  if(~isempty(image_prefixes))
    fprintf(fout_info, 'DIR %d %s\n', case_id, image_sub_dirs{1});
  end
  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))
      continue;
    end
    image_prefix = image_prefixes{tile_id};
    
    if(config.is_verbose)
      fprintf('tile_id: %s\n', image_prefixes{tile_id});
    end
    if(~isfield(config.align.global_align, 'storage_format') || ...
        strcmp(config.align.global_align.storage_format, 'mat')==1)
      load2([tform_dir, image_prefixes{tile_id}, ...
        '.patchwise_affine_transforms.mat'], 'transform');
    else
      fin_t = fopen([tform_dir, image_prefixes{tile_id}, ...
        '.global_patchwise_tform.txt'], 'rt');
      transform.transform = fscanf(fin_t, '%g', [6 inf]);
      fclose(fin_t);
    end      

    image_name = [get_stack_dir(config), image_prefix, '.tif'];
    
    % fold mask file name
    file_name_prefix = [get_fold_dir(config), image_prefix, '.fold_mask'];
    if(isfield(config.fold, 'save_as_tif') && ~isempty(config.fold.save_as_tif) ...
        && config.fold.save_as_tif)
      if(isfield(config.fold, 'save_suffix'))
        file_name_prefix = strrep(file_name_prefix, '.fold_mask', config.fold.save_suffix);
      end
    else
      error('Meant for fold masks saved as .tif files');
    end
    file_name = [file_name_prefix, '.tif'];
    fprintf(fout_info, 'FOLDMAP %s %s\n', image_name, file_name);
    
    % superpixel map
    if(isfield(config, 'segmentation_2D'))
      config_segmentation_2D = config.segmentation_2D;
      is_replaced = true;
    else
      is_replaced = false;
    end
    config.segmentation_2D = config.superpixel_2_seg(1);
    [sp_method, sp_suffixes] = get_superpixel_suffixes(config, image_prefixes{tile_id});
    sp_dir = [get_reconstruction_dir(config), config.segmentation_2D.dir, sp_method, '/'];
    if(is_replaced)
       config.segmentation_2D = config_segmentation_2D;
    end
    file_name = get_storage_file_name([sp_dir, image_prefixes{tile_id}, ...
      sp_suffixes{1},'.mat']);
    sp = load(file_name);
    imwrite(uint16(sp.label_map), [file_name(1:end-4), '.png'], 'BitDepth', 16);
    fprintf(fout_info, 'SPMAP %s %s\n', image_name, [file_name(1:end-4), '.png']);
    
    % copy the transform to the data_to_be_proofread dir for future
    % record.
    check_for_dir([tform_save_dir, image_sub_dirs{tile_id}]);
    for patch_id = 1:size(transform.transform,2)
      fprintf(fout_info, 'TRANSFORM %s::%d %g %g %g %g %g %g\n', ...
        image_name, patch_id, transform.transform(:, patch_id));
    end
      
    % boundary map file name
    file_name_prefix = [get_filtered_image_dir(config), image_prefix];
    file_name = [file_name_prefix, config.split3D.filter_version, '.tif'];
    fprintf(fout_info, 'BOUNDARYMAP %s %s\n', image_name, file_name);
  end
end
fclose(fout_info);

if(config.is_verbose)
  fprintf('STOP: generate_matlab_global_tforms_multi_tile_patches\n');
end
return
end
