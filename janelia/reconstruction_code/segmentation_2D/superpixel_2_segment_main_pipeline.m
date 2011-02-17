function superpixel_2_segment_main_pipeline(config, main_module_id, sub_module_id)
if(config.is_verbose)
  fprintf('START: superpixel_2_segment_main_pipeline\n');
end
switch(main_module_id)
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Superpixel to segment
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 550
    superpixel_2_segment_pipeline(config, sub_module_id)

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interactively choose segmentation parameters for each image
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 590
    gui_choose_segmentation_param(config);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set segmentation parameters for each image as specified in
    % config.segmentation_choose.set_choice
    % This is useful if blocks of sections have different parameters.
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 591
    set_segmentation_suffix_choice(config);
end
if(config.is_verbose)
  fprintf('STOP: superpixel_2_segment_main_pipeline\n');
end
return
end
