function preprocess_images_pipeline(config, main_module_id)

if(config.is_verbose)
  fprintf('START: preprocess_images_pipeline\n');
end
if(~isfield(config, 'precompute') || isempty(config.precompute))
  if(config.is_verbose)
    fprintf('STOP: preprocess_images_pipeline\n');
  end
  return
end

switch(main_module_id)
  case 10
    compute_save_fold_masks(config);
    
  case 20
    dump_BEL_training_data(config); % TODO
  case 21
    boundary_BEL_dump_testing_data_SAB(config);
  case 22
    boundary_BEL_copy_test_results_to_boundary_maps(config);
   
  case 30
    compute_filtered_images(config);
end

if(config.is_verbose)
  fprintf('STOP: preprocess_images_pipeline\n');
end
return
end
