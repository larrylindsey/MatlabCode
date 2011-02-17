function superpixel_2_segment_pipeline(config, sub_module_id)
switch(sub_module_id)
  case 0
    superpixel_config = config.superpixel;
    config.superpixel = superpixel_config(end);
    replace_flag = 0;
    if(isfield(config, 'segmentation_2D'))
      segmentation_config = config.segmentation_2D;
      replace_flag = 1;
    end;
    fprintf('--- Superpixel 2 segment stage %d ---\n', 1);
    config.segmentation_2D = config.superpixel_2_seg(1);
    superpixel_2_segment_main(config);
    config = rmfield(config, 'segmentation_2D');
    for i = 2:length(config.superpixel_2_seg)
      fprintf('--- Superpixel 2 segment stage %d ---\n', i);
      config.segmentation_2D = config.superpixel_2_seg(i);
      segment_2_segment_main(config);
      config = rmfield(config, 'segmentation_2D');
    end
    if(replace_flag==1)
      config.segmentation_2D = segmentation_config;
    end;
    config.superpixel = superpixel_config;
  otherwise
    segmentation_stage = floor(sub_module_id/10);
    fprintf('--- Superpixel 2 segment stage %d ---\n', segmentation_stage);
    superpixel_config = config.superpixel;
    config.superpixel = superpixel_config(end);
    replace_flag = 0;
    if(isfield(config, 'segmentation_2D'))
      segmentation_config = config.segmentation_2D;
      replace_flag = 1;
    end;
    config.segmentation_2D = config.superpixel_2_seg(segmentation_stage);
    if(segmentation_stage==1)
      superpixel_2_segment_main(config);
    else
      segment_2_segment_main(config);
    end
    config = rmfield(config, 'segmentation_2D');
    if(replace_flag==1)
      config.segmentation_2D = segmentation_config;
    end;
    config.superpixel = superpixel_config;
end
return
end
