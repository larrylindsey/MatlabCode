function segmentation_stich_pipeline(config, main_module_id, sub_module_id) %#ok<INUSD>
if(config.is_verbose)
  fprintf('START: segmentation_stich_pipeline\n');
end
switch(main_module_id)
  case 650
    replace_flag = 0;
    if(isfield(config, 'segmentation_2D'))
      segmentation_config = config.segmentation_2D;
      replace_flag = 1;
    end;
    config.segmentation_2D = config.superpixel_2_seg(end);
    align_segment_map_multi_tile(config);
    config = rmfield(config, 'segmentation_2D');
    if(replace_flag==1)
      config.segmentation_2D = segmentation_config;
    end;
  case 660
    % generate identity mapping for locked labels.
    % useful when combining proofread volumes.
    generate_identity_align_segment_mappings(config);
end
if(config.is_verbose)
  fprintf('STOP: segmentation_stich_pipeline\n');
end
return
end
