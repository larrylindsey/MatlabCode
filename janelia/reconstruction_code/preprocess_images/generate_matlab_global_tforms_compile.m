function generate_matlab_global_tforms_compile(config)
% generate_matlab_global_tforms(config)
% generate MATLAB style transforms. These transforms
% will be used for rendering al and cat.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
%

if(isfield(config.stack, 'image_structure') && ~isempty(config.stack.image_structure))
  if(~isfield(config.stack, 'fold') || ...
      ~isfield(config.stack.fold, 'is_considered_in_segmentation') || ...
      ~config.stack.fold.is_considered_in_segmentation)
    error('TO DO');
  else
    generate_matlab_global_tforms_multi_tile_patches_compile(config);
  end
else
  fprintf('Uni tile stack, no global transforms to generate.\n');
end

return
end
