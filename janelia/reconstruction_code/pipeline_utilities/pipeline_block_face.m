function pipeline_block_face(config, main_module_id, sub_module_id)
% pipeline_block_face(module_ids, config)
% Calls different modules of the pipeline based on module_id. To be used
% for Focus Ion Beam or Block Face EM Images.
%
% module_ids: 1000-1999.99
%   1000    precomputation *
%   1001    generate image sub-stack TIFFs after filtering
%   1002    mitochondria processing *
%   1100    initial 3D watershed
%   1110    watershed to superpixels
%   1120    superpixel to segments *
%   1700    combine segment maps of sub-stacks. Useful if stack is very
%             large.
%   1900    generate TIFF stacks of color mapped label for viewing *
%   1910    output data for proofreading *
%
% * To be implemented
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04162008  init code
% v1  01192009  pipeline version 3
% v2  03112009  branch into 2D and 3D reconstruction
%

switch(main_module_id)
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Generate image sub-stacks. Useful if the stack is too large
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 1001
    generate_image_stack(config);

  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute initial 3D segmentationm, e.g., using watershed. Superpixels are
  % agglomerations of this segmentation.
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 1100
    replace_flag = false;
    if(isfield(config, 'segmentation_3D'))
      config_segmentation_3D = config.segmentation_3D;
      replace_flag = true;
    end
    config.segmentation_3D = config.initial_segmentation_3D;
    segment_3D_watershed(config);
    if(replace_flag)
      config.segmentation_3D = config_segmentation_3D;
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute superpixels from initial watershed
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 1110
    replace_flag = false;
    if(isfield(config, 'segmentation_3D'))
      config_segmentation_3D = config.segmentation_3D;
      replace_flag = true;
    end
    config.segmentation_3D = config.superpixel_3D;
    superpixel_2_segment_3D_ladder(config);
    if(replace_flag)
      config.segmentation_3D = config_segmentation_3D;
    end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Combine segment maps of sub-stacks. Useful if stack is very large.
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 1700
    combine_segment_maps_3D(config);
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % Generate TIFF stacks of color mapped label for viewing
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  case 1900
    generate_segment_display_stack(config);
    
  otherwise
    error('Pipeline module id not understood.');
end;

return
end
