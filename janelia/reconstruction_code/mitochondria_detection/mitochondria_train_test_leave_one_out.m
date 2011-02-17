function mitochondria_train_test_leave_one_out(config)
% mitochondria_train_test_leave_one_out(config)
% train and test using Leave-one-out protocol.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

config_0 = config;
for t = config_0.mitochondria.train.case_ids
  config.stack.case_ids = t;
  config.stack.roi = [];
  config.mitochondria.train.case_ids = setdiff(config_0.mitochondria.train.case_ids, t);
  config.stack.image_prefix = config_0.mitochondria.train.image_prefix;
  pipeline([3 4], config);
end

return
end
