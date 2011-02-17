function mitochondria_segment_label_threshold(config)
% mitochondria_segment_label_threshold(config)
% Threshold mitchondria confidence values
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

if(~isfield(mitochondria_config, 'is_verbose'))
  mitochondria_config.is_verbose = true;
end

mito_dir = [get_reconstruction_dir(config), mitochondria_config.dir];

if(mitochondria_config.is_verbose)
  fprintf('\nSTART: mitochondria_segment_label_threshold\n');
end

for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  if(mitochondria_config.is_verbose)
    fprintf('plane:%d\n', case_id);
  end
  image_prefixes = get_image_prefixes_subdirs(config, case_id);
  image_prefix = image_prefixes{1};
if(mitochondria_config.is_verbose)
    fprintf('Loading segment confidence values ... ');
end

%%%% load segment confidence values for image
load_file_name = [mito_dir, image_prefix,...
  '.mito.seg', segment_config.mito_seg_suffix, ...
  '.features', feature_config.features_suffix, ...
  '.gt', segment_config.gt_suffix, ...
  '.model', model_config.model_suffix, ...
  '.conf_values',...
  '.mat'];
seg_conf = load2(load_file_name, 'seg_mt_conf_vals');
num_segs=length(seg_conf.seg_mt_conf_vals);

if(mitochondria_config.is_verbose)
    fprintf('done.\n');
end
 
if(mitochondria_config.is_verbose)
    fprintf('Applying Threshold ... ');
end

%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
%seg_conf_thresh=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seg_labels = (seg_conf.seg_mt_conf_vals>classify_config.seg_conf_thresh);

  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end

  if(mitochondria_config.is_verbose)
    fprintf('Saving mitochondria segments ... ');
  end
  
  %% Save mt segments
  save_file_name = [mito_dir, image_prefix,...
    '.mito.seg', segment_config.mito_seg_suffix, ...
    '.features', feature_config.features_suffix, ...
    '.gt', segment_config.gt_suffix, ...
    '.model', model_config.model_suffix, ...
    '.conf_threshold.', num2str(classify_config.seg_conf_thresh), ...
    '.mat'];
  save2(save_file_name, 'seg_labels'); 
  
end
  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end
  
if(mitochondria_config.is_verbose)
  fprintf('\nSTOP: mitochondria_segment_label_threshold\n');
end
return;
end
