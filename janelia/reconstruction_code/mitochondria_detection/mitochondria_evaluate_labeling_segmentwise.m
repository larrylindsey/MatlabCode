function mitochondria_evaluate_labeling_segmentwise(config)
% mitochondria_evaluate_labeling_segmentwise(config)
% Evaluates detections segmentwise
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
  fprintf('\nSTART: mitochondria_evaluate_labeling_segmentwise\n');
end

results.precision=zeros(length(stack_config.case_ids),1);
results.recall=zeros(length(stack_config.case_ids),1);

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
  segment_config.evaluate_labeling_suffix, ...
  '.mat'];
mt_labels = load2(load_file_name, 'seg_labels');

  load_file_name = [mito_dir, image_prefix, ...
    '.mito.seg', segment_config.mito_seg_suffix, ...
    '.gt', segment_config.gt_suffix, ...
    '.mat'];
  gt_labels = load2(load_file_name, 'seg_labels');
  
  load_file_name = [mito_dir, image_prefix, ...
    '.mito.seg', segment_config.mito_seg_suffix, ...
    '.features', feature_config.features_suffix, ...
    '.mat'];
  segment_features = load2(load_file_name, 'features');
  areas=segment_features.features(:,1);
  
  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end
  
   if(mitochondria_config.is_verbose)
    fprintf('Computing PR values ... ');
    end
  
  num_pixels_correct = sum(mt_labels.seg_labels.*gt_labels.seg_labels.*areas);
  num_pixels_detected = sum(mt_labels.seg_labels.*areas);
  num_true_pixels = sum(gt_labels.seg_labels.*areas);
  
  results.precision(i) = num_pixels_correct / num_pixels_detected * 100;
  results.recall(i) = num_pixels_correct / num_true_pixels * 100;
  
save_file_name = [mito_dir, image_prefix,...
  segment_config.evaluate_labeling_suffix, ...
  '.PRresults.mat'];
save2(save_file_name, 'results');
  
    if(mitochondria_config.is_verbose)
    fprintf('done.\n');
    end
end
  results.precision
  results.recall
  mean(results.precision)
  mean(results.recall)
  if(mitochondria_config.is_verbose)
  fprintf('\nSTOP: mitochondria_evaluate_labeling_segmentwise\n');
end
return;
end
  
  