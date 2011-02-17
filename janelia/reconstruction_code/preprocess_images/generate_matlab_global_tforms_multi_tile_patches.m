function generate_matlab_global_tforms_multi_tile_patches(config)
% generate_matlab_global_tforms_multi_tile_patches(config)
% generate MATLAB style transforms for all the patches when there are
% multiple tiles and patches per section (trakEM xml). These transforms
% will be used for rendering al and cat.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
% v1  09302008  modified to multiple tile per section (trakEM xml)
%

if(config.is_verbose)
  fprintf('START: generate_matlab_global_tforms_multi_tile_patches\n');
end

stack_config = config.stack;
switch(config.align.global_align.method)
  case 'SIFT'
    tform_dir = [get_region_dir(config), config.SIFT.dir];
    tform_save_dir = [get_to_be_proofread_dir(config), config.SIFT.dir];
  case 'deformable_mesh'
    tform_dir = [get_region_dir(config), config.deformable_mesh.dir];
    tform_save_dir = [get_to_be_proofread_dir(config), ...
      config.deformable_mesh.dir];
end

check_for_dir(get_global_stitching_param_dir(config));

images = get_image_from_stack(config, stack_config.case_ids(1));
[height, width] = size(images{1});
for layer_id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(layer_id);
  fprintf('plane: %d\n', case_id);
  
  tforms = {};
  xdata = {};
  ydata = {};
  
  [image_prefixes, image_sub_dirs]  = ...
    get_image_prefixes_subdirs(config, case_id);

  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))
      continue;
    end
    if(config.is_verbose)
      fprintf('tile: %s\n', image_prefixes{tile_id});
    end
    if(~isfield(config.align.global_align, 'storage_format') || ...
        strcmp(config.align.global_align.storage_format, 'mat')==1)
      load2([tform_dir, image_prefixes{tile_id}, ...
        '.patchwise_affine_transforms.mat'], 'transform');
      
      % copy the transform to the data_to_be_proofread dir for future
      % record.
      check_for_dir([tform_save_dir, image_sub_dirs{tile_id}]);
      save2([tform_save_dir, image_prefixes{tile_id}, ...
        '.patchwise_affine_transforms.mat'], 'transform');
      
      n_patch = size(transform.transform, 2);
      for patch_id = 1:n_patch
        if(config.is_verbose)
          fprintf('patch: %d\n', patch_id);
        end
        transform_affine = reshape(transform.transform(:,patch_id), [2 3]);
        tforms{tile_id, patch_id} = ...
          maketform('affine', transform_affine'); %#ok<AGROW>
        [junk, xdata{tile_id, patch_id}, ydata{tile_id, patch_id}] = ...
          imtransform(ones(height, width), tforms{tile_id, patch_id}, 'nearest'); %#ok<ASGLU,AGROW>
      end
    else
      fin_t = fopen([tform_dir, image_prefixes{tile_id}, ...
        '.global_patchwise_tform.txt'], 'rt');
      transform.transform = fscanf(fin_t, '%g', [6 inf]);
      fclose(fin_t);
      
      % copy the transform to the data_to_be_proofread dir for future
      % record.
      check_for_dir([tform_save_dir, image_sub_dirs{tile_id}]);
      fout_t = fopen([tform_save_dir, image_prefixes{tile_id}, ...
        '.global_patchwise_tform.txt'], 'wt');
      fprintf(fout_t, '%g %g %g %g %g %g\n', transform.transform);
      fclose(fout_t);
      
      n_patch = size(transform.transform, 2);
      for patch_id = 1:n_patch
        if(config.is_verbose)
          fprintf('patch: %d\n', patch_id);
        end
        transform_affine = reshape(transform.transform(:,patch_id), [3 2])';
        tforms{tile_id, patch_id} = ...
          maketform('affine', transform_affine'); %#ok<AGROW>
        [junk, xdata{tile_id, patch_id}, ydata{tile_id, patch_id}] = ...
          imtransform(ones(height, width), tforms{tile_id, patch_id}, 'nearest'); %#ok<ASGLU,AGROW>
      end
    end
    
    % copy fold mask to data-to-be-proofread transforms directory
    file_name_prefix = [get_fold_dir(config), image_prefixes{tile_id}, '.fold_mask'];
    fold_mask = load_fold_mask(file_name_prefix, config);
    imwrite(fold_mask, [tform_save_dir, image_prefixes{tile_id}, '.fold_mask.tif']);
  end
  
  save2([get_global_stitching_param_dir(config), ...
    'global_stitching_parameters.', num2str(case_id), ...
    '.mat'], 'xdata', 'ydata', 'tforms');
end

if(config.is_verbose)
  fprintf('STOP: generate_matlab_global_tforms_multi_tile_patches\n');
end
return
end
