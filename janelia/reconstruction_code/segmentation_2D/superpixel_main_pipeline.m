function superpixel_main_pipeline(config, main_module_id, sub_module_id)
if(config.is_verbose)
  fprintf('START: superpixel_main_pipeline\n');
end
switch(main_module_id)
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Train segmentation algorithms
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 410
    switch(sub_module_id)
      case 1
        collect_segment_boundary_sets_groundtruth(config);
      case 2
        collect_segment_boundary_sets_auto_segment_grayscale(config);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Superpixel
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 450
    superpixel_pipeline(config, sub_module_id);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interactively choose superpixel parameters for each image
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 490
    replace_flag = 0;
    if(isfield(config, 'segmentation_choose'))
      segmentation_choose_config = config.segmentation_choose;
      replace_flag = 1;
    end;
    config.segmentation_choose = config.superpixel_choose;
    gui_choose_segmentation_param(config);
    if(replace_flag==1)
      config.segmentation_choose = segmentation_choose_config;
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set superpixel parameters for each image as specified in
    % config.superpixel_choose.set_choice
    % This is useful if blocks of sections have different parameters.
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 491
    set_superpixel_suffix_choice(config);
end
if(config.is_verbose)
  fprintf('STOP: superpixel_main_pipeline\n');
end
return
end
