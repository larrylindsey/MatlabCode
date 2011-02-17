function compile_linkage_evaluate_segment_rand_error(config)
% compile_linkage_evaluate_segment_rand_error(config)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%

if(isfield(config.stack, 'image_structure') && ...
    ~isempty(config.stack.image_structure))
  compile_linkage_evaluate_segment_rand_error_multi_tile(config);
else
  error('No need for compilation');
end

return
end
