function generate_al(config)
% generate_al(config)
% generate al
%
% al{z}'s are the EM images used for the reconstruction
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
%

if(config.is_verbose)
  fprintf('START: generate_al\n');
end
if(isfield(config.stack, 'image_structure') && ~isempty(config.stack.image_structure))
  if(~isfield(config.stack, 'fold') || ...
      ~isfield(config.stack.fold, 'is_considered_in_segmentation') || ...
      ~config.stack.fold.is_considered_in_segmentation)
    generate_al_multi_tile(config);
  else
    generate_al_multi_tile_patches(config);
  end
else
  generate_al_uni_tile(config);
end
if(config.is_verbose)
  fprintf('STOP: generate_al\n');
end
return
end
