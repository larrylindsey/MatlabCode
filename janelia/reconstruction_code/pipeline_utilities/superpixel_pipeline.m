function superpixel_pipeline(config, sub_module_id)
if(config.is_verbose)
  fprintf('START: superpixel_pipeline\n');
end
switch(sub_module_id)
  case 0
    fprintf('--- Superpixel segmentation stage %d ---\n', 1);
    replace_flag = 0;
    if(isfield(config, 'segmentation_2D'))
      segmentation_config = config.segmentation_2D;
      replace_flag = 1;
    end;
    config.segmentation_2D = config.superpixel(1);
    segment_2D_main(config);
    config = rmfield(config, 'segmentation_2D');
    superpixel_config = config.superpixel;
    for i = 2:length(config.superpixel)
      fprintf('--- Superpixel segmentation stage %d ---\n', i);
      config.superpixel = superpixel_config(i-1);
      config.segmentation_2D = superpixel_config(i);
      superpixel_2_superpixel_main(config);
      config = rmfield(config, 'segmentation_2D');
    end
    if(replace_flag==1)
      config.segmentation_2D = segmentation_config;
    end;
    config.superpixel = superpixel_config;
  otherwise
    segmentation_stage = floor(sub_module_id/10);
    fprintf('--- Superpixel segmentation stage %d ---\n', segmentation_stage);
    replace_flag = 0;
    if(isfield(config, 'segmentation_2D'))
      segmentation_config = config.segmentation_2D;
      replace_flag = 1;
    end;
    superpixel_config = config.superpixel;
    if(segmentation_stage==1)
      config.segmentation_2D = superpixel_config(segmentation_stage);
      segment_2D_main(config);
    else
      config.superpixel = superpixel_config(segmentation_stage-1);
      config.segmentation_2D = superpixel_config(segmentation_stage);
      superpixel_2_superpixel_main(config);
    end
    config = rmfield(config, 'segmentation_2D');
    if(replace_flag==1)
      config.segmentation_2D = segmentation_config;
    end;
    config.superpixel = superpixel_config;
end
if(config.is_verbose)
  fprintf('STOP: superpixel_pipeline\n');
end
return
end
