function mitochondria_display_3D_volume(config)
% mitochondria_display_3D_volume(config)
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
  fprintf('\nSTART: mitochondria_display_3D_volume\n');
end


mt_vol=zeros(863,863,length(stack_config.case_ids));

if(segment_config.display_gt)
  gt_vol=zeros(863,863,length(stack_config.case_ids));
end

for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  if(mitochondria_config.is_verbose)
    fprintf('plane:%d\n', case_id);
  end
  if(mitochondria_config.is_verbose)
    fprintf('Loading mitochondria detections ... ');
  end

  image_prefixes = get_image_prefixes_subdirs(config, case_id);
  image_prefix = image_prefixes{1};
  load_file_name = [mito_dir, image_prefix,...
    segment_config.display_suffix, ...
    '.mat'];
  mitochondria = load2(load_file_name, 'output');
  mt_vol(:,:,i) = mitochondria.output;

  if(segment_config.display_gt)
    ground_truth = load2([mito_dir, mitochondria_config.train.dir, image_prefix, ...
      '.mitochondria_annot.mat'], 'zone_mask');
    ground_truth.zone_mask=logical(ground_truth.zone_mask);

    ground_truth.zone_mask(:,1:mitochondria_config.border-1)=0;
    ground_truth.zone_mask(:,end-mitochondria_config.border:end)=0;
    ground_truth.zone_mask(1:mitochondria_config.border-1,:)=0;
    ground_truth.zone_mask(end-mitochondria_config.border:end,:)=0;

    gt_vol(:,:,i) = ground_truth.zone_mask;
  end

end

if(mitochondria_config.is_verbose)
  fprintf('done.\n');
end

if(mitochondria_config.is_verbose)
  fprintf('Making 3D volume ... ');
end

display_volume = reducevolume(mt_vol, [3 3 1]);
figure(9999);
clf;
patches = patch(isosurface(display_volume, 0));
reducepatch(patches, 0.2);
set(patches, 'FaceColor', [1.000, 0.891, 0.009], 'FaceAlpha', 1.0, 'LineStyle', 'none');
hLight = light();

if(segment_config.display_gt)
  display_volume = reducevolume(gt_vol, [3 3 1]);
  figure(8888);
  clf;
  patches = patch(isosurface(display_volume, 0));
  reducepatch(patches, 0.2);
  set(patches, 'FaceColor', [1.000, 0.891, 0.009], 'FaceAlpha', 1.0, 'LineStyle', 'none');
  hLight = light();
end

if(mitochondria_config.is_verbose)
  fprintf('done.\n');
end

if(mitochondria_config.is_verbose)
  fprintf('\nSTOP: mitochondria_display_3D_volume\n');
end
return;
end
