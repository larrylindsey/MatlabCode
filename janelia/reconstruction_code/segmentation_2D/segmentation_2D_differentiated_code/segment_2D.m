function segment_2D(config)
% segment_2D(config)
% Calls a 2D segmentation algorithm depending upon
% config.segmentation_2D.method.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v1  04302008  init code
%

segment_2D_main(config);
return;

switch(config.segmentation_2D.method)
  case 'grayscale_ladder'
    segment_2D_main(config);
  case 'pb_ladder'
    segment_2D_pb_ladder(config);
  case 'BEL_ladder'
    segment_2D_BEL_ladder(config);
  case 'grayscale_AglMeanB'
    segment_2D_grayscale_agglo_mean_boundary(config);
  otherwise
    error('Segmentation algorithm type not recognized');
end

end