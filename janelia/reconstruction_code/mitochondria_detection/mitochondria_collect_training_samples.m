function mitochondria_collect_training_samples(config)
% mitochondria_collect_training_samples(config)
% Calls a function for collecting samples depending upon
% config.mitochondria.feature.type.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v1  05012008  init code
%

switch(config.mitochondria.feature.type)
  case {'heq_intensity_hist', 'intensity_hist'}
    mitochondria_collect_training_samples_intensity_histograms(config);
  case {'heq_intensity_hist_gabor', 'intensity_hist_gabor', ...
      'heq_intensity_hist_heq_gabor'}
    mitochondria_collect_training_samples_intensity_hist_gabor(config);
  case {'heq_intensity_hist_3sec', 'intensity_hist_3sec'}
    mitochondria_collect_training_samples_intensity_histograms_3sec(config);
  case {'heq_intensity_hist_3s2', 'intensity_hist_3s2'}
    mitochondria_collect_training_samples_intensity_histograms_3s2(config);
  case {'heq_intensity_norm_hist_3s2', 'intensity_norm_hist_3s2'}
    mitochondria_collect_training_samples_intensity_norm_hist_3s2(config);
  otherwise
    error('Mitochondria feature type not recognized');
end

end