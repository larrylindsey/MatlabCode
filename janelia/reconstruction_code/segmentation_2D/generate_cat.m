function generate_cat(config)
% generate_cat(config)
% generate cat
%
% cat{z}'s are the superpixel maps being used in the reconstruction
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
%

if(isfield(config.stack, 'image_structure') && ~isempty(config.stack.image_structure))
  generate_cat_multi_tile(config);
else
  generate_cat_uni_tile(config);
end

return
end
