function linkage_pipeline(config, main_module_id, sub_module_id)
if(config.is_verbose)
  fprintf('START: linkage_pipeline\n');
end
switch(main_module_id)
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Linkage - training
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 710
    linkage_3D_train(config, sub_module_id);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Linkage - computation on stack
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 750
    replace_flag = 0;
    if(isfield(config, 'segmentation_2D'))
      segmentation_config = config.segmentation_2D;
      replace_flag = 1;
    end;
    config.segmentation_2D = config.superpixel_2_seg(end);
    linkage_3D_generate_linkage_graph(config);
    config = rmfield(config, 'segmentation_2D');
    if(replace_flag==1)
      config.segmentation_2D = segmentation_config;
    end;
end
if(config.is_verbose)
  fprintf('STOP: linkage_pipeline\n');
end
return
end
