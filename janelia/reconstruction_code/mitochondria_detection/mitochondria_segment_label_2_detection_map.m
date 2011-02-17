function mitochondria_segment_label_2_detection_map(config)
% mitochondria_segment_label_2_detection_map(config)
% Converts segment labels to a detection map
%
% Nicholas Sofroniew,
% Univ. of Cambridge, UK.
% Visitor March 2009, JFRC, HHMI.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  03132009  init code
%

stack_config = config.stack;
mitochondria_config = config.mitochondria;
segment_config = mitochondria_config.segment;
classify_config = segment_config.classify;
feature_config = classify_config.feature;
model_config = classify_config.model;
mrf_config = segment_config.mrf;

if(~isfield(mitochondria_config, 'is_verbose'))
  mitochondria_config.is_verbose = true;
end

mito_dir = [get_reconstruction_dir(config), mitochondria_config.dir];

if(mitochondria_config.is_verbose)
  fprintf('\nSTART: mitochondria_segment_label_2_detection_map\n');
end

for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  if(mitochondria_config.is_verbose)
    fprintf('plane:%d\n', case_id);
  end
  image_prefixes = get_image_prefixes_subdirs(config, case_id);
  image_prefix = image_prefixes{1};

  if(mitochondria_config.is_verbose)
    fprintf('Loading segment labels ... ');
  end

  load_file_name = [mito_dir, image_prefix,...
    segment_config.label2detection_suffix, ...
    '.mat'];
  labels = load2(load_file_name, 'seg_labels');


  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end

  if(mitochondria_config.is_verbose)
    fprintf('Saving segment detections ... ');
  end

  %%% Find segments labeled as mitochondria
  load_file_name = [mito_dir, image_prefix, ...
    '.mito.seg', segment_config.mito_seg_suffix,...
    '.mat'];
  segments = load2(load_file_name, 'label_map');
  mt_detections=ismember(segments.label_map,find(labels.seg_labels==1));

  if(mitochondria_config.is_verbose_figures)
    figure(2*i);
    imshow(segments.label_map);
    title('Mitochondria segmentation');
  end
  if(mitochondria_config.is_verbose_figures)
    figure(2*i+1);
    imshow(mt_detections);
    title('Mitochondria detections segmentwise');
  end

  %% Save mt segments
  save_file_name = [mito_dir, image_prefix,...
    segment_config.label2detection_suffix, ...
    '.detections.mat'];
  save2(save_file_name, 'mt_detections');


  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end

end

if(mitochondria_config.is_verbose)
  fprintf('\nSTOP: mitochondria_segment_label_2_detection_map\n');
end
return;
end