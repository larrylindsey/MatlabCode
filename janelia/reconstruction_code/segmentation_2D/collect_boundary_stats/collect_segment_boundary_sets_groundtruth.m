function collect_segment_boundary_sets_groundtruth(config)
% collect_segment_boundary_sets_groundtruth(config)
% Collect sets of segment boundary elements from groundtruth
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  08272008  init code
%

switch(config.train_segmentation.boundary_method)
  case 'grayscale'
    collect_segment_boundary_sets_groundtruth_grayscale(config);
  otherwise
    error('Boundary method not recognized for collecting boundary sets');
end

end
