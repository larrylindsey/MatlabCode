function linkage_evaluate_segment_rand_error(config)
% linkage_evaluate_segment_rand_error(config)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%

if(isfield(config.stack, 'image_structure') && ...
    ~isempty(config.stack.image_structure))
  linkage_evaluate_segment_rand_error_multi_tile(config);
else
  linkage_evaluate_segment_rand_error_uni_tile(config);
end

return
end
