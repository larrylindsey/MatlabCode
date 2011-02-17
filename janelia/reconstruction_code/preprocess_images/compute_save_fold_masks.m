function compute_save_fold_masks(config)
% compute_save_fold_masks(config)
% Compute fold masks for images
%
%
% Input:
%   config    config datastructure of the reconstruction
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  01022009  init. code
%

fold_config = config.precompute.fold;
if(~isfield(fold_config, 'is_verbose'))
  fold_config.is_verbose = true;
end
if(fold_config.is_verbose)
  fprintf('START: compute_save_fold_masks\n');
end

stack_config = config.stack;
fold_dir = get_fold_dir(config);

for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  if(fold_config.is_verbose)
    fprintf('case_id: %d\n', case_id);
  end
  [images, image_prefixes, image_sub_dirs, is_to_be_processed] = ...
    get_image_from_stack(config, case_id);
  
  for j = 1:length(images)
    image = images{j};
    image_prefix = image_prefixes{j};
    image_sub_dir = image_sub_dirs{j};
    if(fold_config.is_verbose)
      fprintf('%s\n', image_prefix);
    end
    if(is_to_be_processed(j)==0)
      if(fold_config.is_verbose)
        fprintf('is_to_be_processed=false, skipping\n');
      end
      continue;
    end
    
    if(isfield(fold_config, 'filter_version') && ~isempty(fold_config.filter_version))
      image = filter_image(image, fold_config.filter_version);
    end
    
    fold_mask = uint8(compute_fold_mask(uint8(255*image), fold_config));
    
    check_for_dir([fold_dir, image_sub_dir]);
    if(~isfield(config.fold, 'save_as_tif') || ~config.fold.save_as_tif)
      file_name = [fold_dir, image_prefix, '.fold_mask.mat'];
      save2(file_name, 'fold_mask');
    else
      imwrite(fold_mask, [fold_dir, image_prefix, '.fold_mask.tif']);
    end
    
    if(isfield(config.precompute.fold, 'save_patch_tifs') && ...
        config.precompute.fold.save_patch_tifs)
      n_patch = max(fold_mask(:));
      for patch_id = 1:n_patch
        patch_mask = 255*(fold_mask == patch_id);
        file_name = [fold_dir, image_prefix, '.patch_mask.', num2str(patch_id), '.tif'];
        imwrite(patch_mask, file_name);
      end
    end
  end
end

if(fold_config.is_verbose)
  fprintf('STOP: compute_save_fold_masks\n');
end

return
end
