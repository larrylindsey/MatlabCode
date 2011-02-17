function generate_annotations(config)
% generate_annotations(config)
% generate annotations
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

if(config.is_verbose)
  fprintf('START: generate_annotations\n');
end
if(isfield(config.stack, 'image_structure') && ~isempty(config.stack.image_structure))
  if(~isfield(config.stack, 'fold') || ...
      ~isfield(config.stack.fold, 'is_considered_in_segmentation') || ...
      ~config.stack.fold.is_considered_in_segmentation)
    error('To Do');
  else
    generate_annotations_multi_tile_patches(config);
  end
else
  error('To Do');
end
if(config.is_verbose)
  fprintf('STOP: generate_annotations\n');
end
return
end
