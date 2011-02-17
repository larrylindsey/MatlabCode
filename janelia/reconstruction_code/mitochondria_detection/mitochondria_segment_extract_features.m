function mitochondria_segment_extract_features(config)
% mitochondria_segment_extract_features(config)
% Extract features to represent mitochondria segments
%
% Nicholas Sofroniew,
% Univ. of Cambridge, UK.
% Visitor March 2009, JFRC, HHMI.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  03122009  init code
%

stack_config = config.stack;
mitochondria_config = config.mitochondria;
segment_config = mitochondria_config.segment;
classify_config = segment_config.classify;
feature_config = classify_config.feature;

if(~isfield(mitochondria_config, 'is_verbose'))
  mitochondria_config.is_verbose = true;
end
if(strcmp(feature_config.type, 'a_cnch_cnch3_va_vna_va3_vna3')~=1)
  error('Feature type is not a_cnch_cnch3_va_vna_va3_vna3');
end
mito_dir = [get_reconstruction_dir(config), mitochondria_config.dir];
vesicle_dir = [get_reconstruction_dir(config), config.vesicle.dir];
if(mitochondria_config.is_verbose)
  fprintf('\nSTART: mitochondria_segment_extract_features\n');
end
for i = 1:length(stack_config.case_ids)
  case_id_0 = stack_config.case_ids(i);
  if(mitochondria_config.is_verbose)
    fprintf('plane:%d\n', case_id_0);
  end
  if(mitochondria_config.is_verbose)
    fprintf('Loading mitochondria detection maps and vesicle maps ... ');
  end
  image_prefixes = get_image_prefixes_subdirs(config, case_id_0);
  image_prefix_0 = image_prefixes{1};
  mitochondria_conf_0 = load2([mito_dir, image_prefix_0, ...
    '.mitochondria_det_conf', feature_config.mitochondria_confidence_suffix, ...
    '.mat'], 'mitochondria_detection_confidence');
  vesicle_mask_0 = load2([vesicle_dir, image_prefix_0, ...
    '.vesicle_mask', feature_config.vesicle_suffix, '.mat'], ...
    'vesicle_mask');
  if(i<length(stack_config.case_ids))
    case_id_u = stack_config.case_ids(i+1);
    image_prefixes = get_image_prefixes_subdirs(config, case_id_u);
    image_prefix_u = image_prefixes{1};
    mitochondria_conf_u = load2([mito_dir, image_prefix_u, ...
      '.mitochondria_det_conf', feature_config.mitochondria_confidence_suffix, ...
      '.mat'], 'mitochondria_detection_confidence');
    vesicle_mask_u = load2([vesicle_dir, image_prefix_u, ...
      '.vesicle_mask', feature_config.vesicle_suffix, '.mat'], ...
      'vesicle_mask');
  end
  if(i>1)
    case_id_d = stack_config.case_ids(i-1);
    image_prefixes = get_image_prefixes_subdirs(config, case_id_d);
    image_prefix_d = image_prefixes{1};
    mitochondria_conf_d = load2([mito_dir, image_prefix_d, ...
      '.mitochondria_det_conf', feature_config.mitochondria_confidence_suffix, ...
      '.mat'], 'mitochondria_detection_confidence');
    vesicle_mask_d = load2([vesicle_dir, image_prefix_d, ...
      '.vesicle_mask', feature_config.vesicle_suffix, '.mat'], ...
      'vesicle_mask');
  end
  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end

  load_file_name = [mito_dir, image_prefix_0, ...
    '.mito.seg', segment_config.mito_seg_suffix,...
    '.mat'];
  segments = load2(load_file_name, 'label_map');

  if(mitochondria_config.is_verbose)
    fprintf('Computing features ... ');
  end
  numSegs=max(segments.label_map(:));
  features=zeros(numSegs, ...
    1 + 1 + 1 + 1 + 1 + 2*length(feature_config.mito_histogram_bins));

  %%% For each segment compute its features
  for indSeg=1:numSegs
    Region=segments.label_map==indSeg;
    AREA=sum(sum(Region));

    %%% Number of vesicles and cummulative hist of Mt conf values for the segment
    VsNUM=sum(sum(Region.*vesicle_mask_0.vesicle_mask));
    MtVALS = mitochondria_conf_0.mitochondria_detection_confidence(Region);
    MtHist=hist(MtVALS, feature_config.mito_histogram_bins)/AREA*100;
    MtCumHistOrig=cumsum(MtHist);

    %%% Number of vesicles and cummulative hist of Mt conf values for a 3D volume
    %%% around the segment
    VsNUM_sum = VsNUM;
    MtVals_3s = MtVALS';
    section_count = 1;
    if(i<length(stack_config.case_ids))
      VsNUMup = sum(sum(Region.*vesicle_mask_u.vesicle_mask));
      VsNUM_sum = VsNUM_sum + VsNUMup;

      MtVALSup = mitochondria_conf_u.mitochondria_detection_confidence(Region);
      MtVals_3s = [MtVals_3s, MtVALSup'];

      section_count = section_count + 1;
    end
    if(i>1)
      VsNUMdown = sum(sum(Region.*vesicle_mask_d.vesicle_mask));
      VsNUM_sum = VsNUM_sum + VsNUMdown;

      MtVALSdown = mitochondria_conf_d.mitochondria_detection_confidence(Region);
      MtVals_3s = [MtVals_3s, MtVALSdown'];

      section_count = section_count + 1;
    end
    VsNUMavg = VsNUM_sum / section_count;

    MtHist = hist(MtVals_3s, feature_config.mito_histogram_bins)/AREA*100/section_count;
    MtCumHist=cumsum(MtHist);

    features(indSeg,:)=[AREA, VsNUM, VsNUM/AREA*100, VsNUMavg, VsNUMavg/AREA*100, ...
      MtCumHist, MtCumHistOrig];
  end
  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end

  %% Write matrix with segment features for each image
  save_file_name = [mito_dir, image_prefix_0, ...
    '.mito.seg', segment_config.mito_seg_suffix, ...
    '.features.', feature_config.type, ...
    '.', feature_config.version, ...
    feature_config.mitochondria_confidence_suffix, ...
    feature_config.vesicle_suffix, ...
    '.mat'];
  save2(save_file_name, 'features');

end;
if(mitochondria_config.is_verbose)
  fprintf('\nSTOP: mitochondria_segment_extract_features\n');
end
return;
end
