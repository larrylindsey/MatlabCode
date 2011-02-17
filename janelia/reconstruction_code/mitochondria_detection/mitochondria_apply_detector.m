function mitochondria_apply_detector(config)
% mitochondria_apply_detector(config)
% Calls a function for apply the mitochondria detector depending upon
% config.mitochondria.feature.type and config.mitochondria.model.type
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v1  05012008  init code
%

switch(config.mitochondria.feature.type)
  case {'heq_intensity_hist', 'intensity_hist'}
    switch(config.mitochondria.model.type)
      case 'boost'
        mitochondria_apply_boosted_detector_intensity_hist(config);
    end
  case {'heq_intensity_hist_gabor', 'intensity_hist_gabor', ...
      'heq_intensity_hist_heq_gabor'}
    switch(config.mitochondria.model.type)
      case 'boost'
        mitochondria_apply_boosted_detector_intensity_hist_gabor(config);
    end
  case {'heq_intensity_hist_3sec', 'intensity_hist_3sec'}
    switch(config.mitochondria.model.type)
      case 'boost'
        mitochondria_apply_boosted_detector_intensity_hist_3sec(config);
    end
  case {'heq_intensity_hist_3s2', 'intensity_hist_3s2'}
    switch(config.mitochondria.model.type)
      case 'boost'
        mitochondria_apply_boosted_detector_intensity_hist_3s2(config);
    end
  case {'heq_intensity_norm_hist_3s2', 'intensity_norm_hist_3s2'}
    switch(config.mitochondria.model.type)
      case 'boost'
        mitochondria_apply_boosted_detector_intensity_norm_hist_3s2(config);
    end
  otherwise
    error('Mitochondria feature type not recognized');
end

end