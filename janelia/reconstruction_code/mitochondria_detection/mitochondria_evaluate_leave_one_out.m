function mitochondria_evaluate_leave_one_out(config)
% mitochondria_evaluate_leave_one_out(config)
% evaluate mitochondria
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

config.stack.case_ids = config.mitochondria.train.case_ids;
config.stack.roi = [];
config.stack.image_prefix = config.mitochondria.train.image_prefix;
pipeline(4.2, config);

return
end
