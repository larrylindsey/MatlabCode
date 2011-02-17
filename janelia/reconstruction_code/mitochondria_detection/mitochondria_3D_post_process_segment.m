function mitochondria_3D_post_process_segment(config)
% mitochondria_segment_extract_features(config)
% Uses morphological filtering to process detections
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
  fprintf('\nSTART: mitochondria_3D_post_process_segment\n');
end

%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
mt_area_threshold=350;
mt_vol_threshold=2000;
StEl=strel('disk',4);
min_3D_height=4;
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
    segment_config.postprocess_suffix, ...
    '.mat'];
  mitochondria_segs = load2(load_file_name, 'mt_detections');
  mitochondria_seg_volumes(:,:,i) = mitochondria_segs.mt_detections;
  %%%% combine nearby segments
  mitochondria = imclose(mitochondria_segs.mt_detections,OutputStEl);

  mitochondria=bwareaopen(mitochondria,mt_area_threshold);
  mitochondria_volumes(:,:,i)=mitochondria;
end

if(mitochondria_config.is_verbose)
  fprintf('done.\n');
end

if(mitochondria_config.is_verbose)
  fprintf('Filtering segments ... ');
end

mitochondria_volumes=bwareaopen(mitochondria_volumes,mt_vol_threshold);

%%% remove detections of small height
mt_2=min(mitochondria_volumes,mitochondria_seg_volumes);
mt_3=mt_2;
indIm=0;
for i = 1:length(stack_config.case_ids)
  mt_2(:,:,i)=imerode(mt_2(:,:,i),StEl);
end

label_mt=bwlabeln(mt_2);
bounding_box=regionprops(label_mt,'BoundingBox');
bounding_box=cat(1,bounding_box.BoundingBox);
remove=find(bounding_box(:,6)<min_3D_height);

mt_2(ismember(label_mt,remove))=0;

for i = 1:length(stack_config.case_ids)
  mt_3(:,:,i)=imreconstruct(mt_2(:,:,i),mt_3(:,:,i));
end

mitochondria_seg_volumes=min(mitochondria_seg_volumes,mt_3);

if(mitochondria_config.is_verbose)
  fprintf('done.\n');
end

if(mitochondria_config.is_verbose)
  fprintf('Saving mitochondria detections ... ');
end

for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  if(mitochondria_config.is_verbose)
    fprintf('plane:%d\n', case_id);
  end
  image_prefixes = get_image_prefixes_subdirs(config, case_id);
  image_prefix = image_prefixes{1};
  mt_detections=mitochondria_seg_volumes(:,:,i);
  mt_detections = logical(mt_detections);

  save_file_name = [mito_dir, image_prefix,...
    segment_config.postprocess_suffix, ...
    '.processed.detections.mat'];
  save2(save_file_name, 'mt_detections');


  load_file_name = [mito_dir, image_prefix, ...
    '.mito.seg', segment_config.mito_seg_suffix,...
    '.mat'];
  segments = load2(load_file_name, 'label_map');

  mt_segments = segments.label_map.*mt_detections;
  indices = unique(mt_segments);
  seg_labels=zeros(max(segments.label_map(:)),1);
  seg_labels(indices(2:end))=1;

  save_file_name = [mito_dir, image_prefix,...
    segment_config.postprocess_suffix, ...
    '.processed.mat'];
  save2(save_file_name, 'seg_labels');

end

if(mitochondria_config.is_verbose)
  fprintf('done.\n');
end

if(mitochondria_config.is_verbose)
  fprintf('\nSTOP: mitochondria_3D_post_process_segment\n');
end
return;
end
