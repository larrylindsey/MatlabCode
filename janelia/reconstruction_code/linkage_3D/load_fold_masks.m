function [fold_masks, n_patch] = load_fold_masks(image_prefixes, config)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

fold_masks = struct('fold_mask', []);
n_patch = zeros(1, length(image_prefixes));
fold_dir = get_fold_dir(config);
for tile = 1:length(image_prefixes)
  file_name_prefix = [fold_dir, image_prefixes{tile}, '.fold_mask'];
  fold_masks(tile).fold_mask = load_fold_mask(file_name_prefix, config);
  n_patch(tile) = max(fold_masks(tile).fold_mask(:));
end
return
end
