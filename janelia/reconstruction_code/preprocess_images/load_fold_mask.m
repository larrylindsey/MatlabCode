function fold_mask = load_fold_mask(file_name_prefix, config)
% fold_mask = load_fold_mask(file_name_prefix, config)
%
% Loads fold masks

if(isfield(config.fold, 'save_as_tif') && ~isempty(config.fold.save_as_tif) ...
    && config.fold.save_as_tif)
  if(isfield(config.fold, 'save_suffix'))
    file_name_prefix = strrep(file_name_prefix, '.fold_mask', config.fold.save_suffix);
  end
  fold_mask = imread([file_name_prefix, '.tif']);
else
  load2([file_name_prefix, '.mat'], 'fold_mask');
end

return
end
