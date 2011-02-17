function superpixel_2_segment(config)
% segment_2D(config)
% Calls a superpixel-to-segment algorithm depending upon
% config.segmentation_2D.method.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v1  04302008  init code
%

switch(config.segmentation_2D.method)
  case 'grayscale_ladder'
    superpixel_2_segment_main(config);
  case 'BEL_ladder'
    superpixel_2_segment_BEL_ladder(config);
  case 'grayscale_AglMeanB'
    superpixel_2_segment_main(config);
  case 'grayscale_AglBIFL'
    superpixel_2_segment_main(config);
  case 'grayscale_AglMedianB'
    superpixel_2_segment_main(config);
  otherwise
    error('Segmentation algorithm type not recognized');
end

end