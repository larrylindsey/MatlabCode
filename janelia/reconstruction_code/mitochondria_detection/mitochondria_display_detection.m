function mitochondria_display_detection(config)
% mitochondria_display_detection(config)
% Display 3D volumes
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
  fprintf('\nSTART: mitochondria_display_detection\n');
end

for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  if(mitochondria_config.is_verbose)
    fprintf('plane:%d\n', case_id);
  end
  if(mitochondria_config.is_verbose)
    fprintf('Loading mitochondria detections ... ');
  end

  [images, image_prefixes] = get_image_from_stack(config, case_id);
  image = images{1};
  image_prefix = image_prefixes{1};
  load_file_name = [mito_dir, image_prefix,...
    segment_config.display_suffix, ...
    '.mat'];
  mitochondria = load2(load_file_name, 'output');

  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end

  if(mitochondria_config.is_verbose)
    fprintf('Displaying detections ... ');
  end

  image = filter_image(image, segment_config.filter_version);
  image(:,1:mitochondria_config.border-1)=0;
  image(:,end-mitochondria_config.border:end)=0;
  image(1:mitochondria_config.border-1,:)=0;
  image(end-mitochondria_config.border:end,:)=0;

  overlay_mt=repmat(im2double(image), [1 1 3]);
  overlay_mt(:,:,3) = mitochondria.output;
  figure(i), imshow(overlay_mt)
  title('Mitochondria detections');

  if(segment_config.display_gt)
    ground_truth = load2([mito_dir, mitochondria_config.train.dir, image_prefix, ...
      '.mitochondria_annot.mat'], 'zone_mask');
    ground_truth.zone_mask=logical(ground_truth.zone_mask);

    ground_truth.zone_mask(:,1:mitochondria_config.border-1)=0;
    ground_truth.zone_mask(:,end-mitochondria_config.border:end)=0;
    ground_truth.zone_mask(1:mitochondria_config.border-1,:)=0;
    ground_truth.zone_mask(end-mitochondria_config.border:end,:)=0;

    overlay_gt=repmat(im2double(image), [1 1 3]);
    overlay_gt(:,:,3) = ground_truth.zone_mask;
    figure(1000+i), imshow(overlay_gt)
    title('Ground truth detections');
  end
  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end
end

if(mitochondria_config.is_verbose)
  fprintf('\nSTOP: mitochondria_display_detection\n');
end
return;
end
