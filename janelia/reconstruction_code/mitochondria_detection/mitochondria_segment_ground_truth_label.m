function mitochondria_segment_ground_truth_label(config)
% mitochondria_segment_gt_labels(config)
% Assign ground truth labels to segments based on perecent area labeled as
% mitochondira by manual annotation
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

if(~isfield(mitochondria_config, 'is_verbose'))
  mitochondria_config.is_verbose = true;
end

mito_dir = [get_reconstruction_dir(config), mitochondria_config.dir];

if(mitochondria_config.is_verbose)
  fprintf('\nSTART: mitochondria_segment_ground_truth_label\n');
end

for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  if(mitochondria_config.is_verbose)
    fprintf('plane:%d\n', case_id);
  end
  image_prefixes = get_image_prefixes_subdirs(config, case_id);
  image_prefix = image_prefixes{1};
  if(mitochondria_config.is_verbose)
    fprintf('Loading segment label map and ground truth map ... ');
  end
  
  load_file_name = [mito_dir, image_prefix, ...
    '.mito.seg', segment_config.mito_seg_suffix,...
    '.mat'];
  segments = load2(load_file_name, 'label_map');
  
  ground_truth = load2([mito_dir, mitochondria_config.train.dir, image_prefix, ...
    '.mitochondria_annot.mat'], 'zone_mask');
  ground_truth.zone_mask=logical(ground_truth.zone_mask);
  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end

  if(mitochondria_config.is_verbose_figures)
    figure(1);
    imshow(segments.label_map);
    title('Segmentation');
  end
  
  if(mitochondria_config.is_verbose_figures)
    figure(2);
    imshow(ground_truth.zone_mask);
    title('Segmentation');
  end
  
  if(mitochondria_config.is_verbose)
    fprintf('Computing segment ground truth labels ... ');
  end

  %%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
  %%  gt_area_threshold=.30;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%% For each segment compute its fraction area overlap with ground truth
  num_segs=max(segments.label_map(:));
  seg_labels=zeros(num_segs,1);
  for seg_index=1:num_segs
    region=segments.label_map==seg_index;
    area=sum(sum(region));
    area_gt=sum(sum(region.*ground_truth.zone_mask));
    fraction_gt=area_gt/area;
    seg_labels(seg_index)=fraction_gt>=segment_config.gt_area_threshold;
  end

  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end

  %% Write matrix with segment ground truth labels for each image
  save_file_name = [mito_dir, image_prefix, ...
    '.mito.seg', segment_config.mito_seg_suffix, ...
    '.gt.', num2str(segment_config.gt_area_threshold), ...
    '.mat'];
  save2(save_file_name, 'seg_labels');
  

end;
if(mitochondria_config.is_verbose)
  fprintf('\nSTOP: mitochondria_segment_ground_truth_label\n');
end
return;
end
