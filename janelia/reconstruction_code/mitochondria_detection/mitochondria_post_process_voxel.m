function mitochondria_post_process_voxel(config)
% mitochondria_post_process_voxel(config)
% Combine nearby segments to get final detections
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
model_config = classify_config.model;
mrf_config = segment_config.mrf;

if(~isfield(mitochondria_config, 'is_verbose'))
  mitochondria_config.is_verbose = true;
end

mito_dir = [get_reconstruction_dir(config), mitochondria_config.dir];

if(mitochondria_config.is_verbose)
  fprintf('\nSTART: mitochondria_post_process_voxel\n');
end

%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
OutputStEl=strel('square',2);
size_thicken=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Should incorporate into config file, filters, and file suffix

mitochondria_volumes=zeros(863,863,length(stack_config.case_ids));
mitochondria_seg_volumes=zeros(863,863,length(stack_config.case_ids));


for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  if(mitochondria_config.is_verbose)
    fprintf('plane:%d\n', case_id);
  end
  if(mitochondria_config.is_verbose)
    fprintf('Loading mitochondria segment detections ... ');
  end

  image_prefixes = get_image_prefixes_subdirs(config, case_id);
  image_prefix = image_prefixes{1};
  load_file_name = [mito_dir, image_prefix,...
    segment_config.voxelpostprocess_suffix, ...
    '.mat'];
  segments = load2(load_file_name, 'mt_detections');

  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end

  if(mitochondria_config.is_verbose)
    fprintf('Process segments ... ');
  end
  %%%% combine nearby segments
  output=imclose(segments.mt_detections,OutputStEl);
  output=bwmorph(output,'thicken',size_thicken);

  save_file_name = [mito_dir, image_prefix,...
    segment_config.voxelpostprocess_suffix, ...
    '.voxelprocessed.mat'];
  save2(save_file_name, 'output');

  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end
end

if(mitochondria_config.is_verbose)
  fprintf('\nSTOP: mitochondria_post_process_voxel\n');
end
return;
end