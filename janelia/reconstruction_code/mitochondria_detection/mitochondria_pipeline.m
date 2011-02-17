function mitochondria_pipeline(config, module_id, sub_module_id)

if(config.is_verbose)
  fprintf('START: mitochondria_pipeline\n');
end
switch(module_id)
  case 210
    mitochondria_manual_annotation(config);
  case 220
    mitochondria_collect_training_samples(config);
  case 230
    mitochondria_train_boosted_detector(config);
  case 236
    mitochondria_train_test_leave_one_out(config);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mitochondria - detection
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 240
    mitochondria_apply_detector(config);
  case 246
    mitochondria_evaluate_detector(config);
  case 247
    mitochondria_evaluate_leave_one_out(config);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mitochondria detection - Nicholas
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 299
    mitochondria_pipeline_NJS(config, sub_module_id);
end
if(config.is_verbose)
  fprintf('STOP: mitochondria_pipeline\n');
end
return
end
