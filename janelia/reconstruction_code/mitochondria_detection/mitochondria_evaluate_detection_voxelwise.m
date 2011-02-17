function mitochondria_evaluate_detection_voxelwise(config)
% mitochondria_evaluate_detection_voxelwise(config)
% Evaluates detections voxelwise
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
  fprintf('\nSTART: mitochondria_evaluate_detection_voxelwise\n');
end

results.precision=zeros(length(stack_config.case_ids),1);
results.recall=zeros(length(stack_config.case_ids),1);
results.precision_5=zeros(length(stack_config.case_ids),1);
results.recall_5=zeros(length(stack_config.case_ids),1);
results.precision_10=zeros(length(stack_config.case_ids),1);
results.recall_10=zeros(length(stack_config.case_ids),1);

for i = 1:length(stack_config.case_ids)
    case_id = stack_config.case_ids(i);
    if(mitochondria_config.is_verbose)
        fprintf('plane:%d\n', case_id);
    end
    image_prefixes = get_image_prefixes_subdirs(config, case_id);
    image_prefix = image_prefixes{1};
    
    if(mitochondria_config.is_verbose)
    fprintf('Loading detections and ground truth ... ');
    end

load_file_name = [mito_dir, image_prefix,...
  segment_config.evaluate_detection_suffix, ...
  '.mat'];
mitochondria = load2(load_file_name, 'output');
  mitochondria.output=logical(mitochondria.output);
  
  ground_truth = load2([mito_dir, mitochondria_config.train.dir, image_prefix, ...
    '.mitochondria_annot.mat'], 'zone_mask');
  ground_truth.zone_mask=logical(ground_truth.zone_mask);
  
   ground_truth.zone_mask(:,1:mitochondria_config.border-1)=0;
    ground_truth.zone_mask(:,end-mitochondria_config.border:end)=0;
    ground_truth.zone_mask(1:mitochondria_config.border-1,:)=0;
    ground_truth.zone_mask(end-mitochondria_config.border:end,:)=0;
  
  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end
  
   if(mitochondria_config.is_verbose)
    fprintf('Computing PR values ... ');
   end
  
   false_accept_map = (~ground_truth.zone_mask) & mitochondria.output;
   false_reject_map = (~mitochondria.output) & ground_truth.zone_mask;
   correct_map = mitochondria.output & ground_truth.zone_mask;
   
  num_pixels_correct = sum(correct_map(:));
  num_pixels_detected = sum(mitochondria.output(:));
  num_true_pixels = sum(ground_truth.zone_mask(:));
  
  results.precision(i) = num_pixels_correct / num_pixels_detected * 100;
  results.recall(i) = num_pixels_correct / num_true_pixels * 100;
  
   accept_distance_map = bwdist(ground_truth.zone_mask);
   reject_distance_map = bwdist(mitochondria.output);
   
   false_accept_distance_map_5=(accept_distance_map<=5).*false_accept_map;
   false_reject_distance_map_5=(reject_distance_map<=5).*false_reject_map;
   
  num_pixels_detected_5 = sum(mitochondria.output(:))-sum(false_accept_distance_map_5(:));
  num_true_pixels_5 = sum(ground_truth.zone_mask(:))-sum(false_reject_distance_map_5(:));
   
  results.precision_5(i) = num_pixels_correct / num_pixels_detected_5 * 100;
  results.recall_5(i) = num_pixels_correct / num_true_pixels_5 * 100;
   
    false_accept_distance_map_10=(accept_distance_map<=10).*false_accept_map;
   false_reject_distance_map_10=(reject_distance_map<=10).*false_reject_map;
   
  num_pixels_detected_10 = sum(mitochondria.output(:))-sum(false_accept_distance_map_10(:));
  num_true_pixels_10 = sum(ground_truth.zone_mask(:))-sum(false_reject_distance_map_10(:));
  
  results.precision_10(i) = num_pixels_correct / num_pixels_detected_10 * 100;
  results.recall_10(i) = num_pixels_correct / num_true_pixels_10 * 100;
  
  
save_file_name = [mito_dir, image_prefix,...
  segment_config.evaluate_detection_suffix, ...
  '.PRvoxelresults.mat'];
save2(save_file_name, 'results');
  
    if(mitochondria_config.is_verbose)
    fprintf('done.\n');
    end
end
  results.precision
  results.recall
  mean(results.precision)
  mean(results.recall)
    results.precision_5
  results.recall_5
  mean(results.precision_5)
  mean(results.recall_5)
    results.precision_10
  results.recall_10
  mean(results.precision_10)
  mean(results.recall_10)
  if(mitochondria_config.is_verbose)
  fprintf('\nSTOP: mitochondria_evaluate_labeling_segmentwise\n');
end
return;
end
  
  