function mitochondria_segment_apply_classifier(config)
% mitochondria_segment_apply_classifier(config)
% Apply boosted classifier to get mitochondria confidence values
% for segments.
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

if(~isfield(mitochondria_config, 'is_verbose'))
  mitochondria_config.is_verbose = true;
end
if(strcmp(feature_config.type, 'a_cnch_cnch3_va_vna_va3_vna3')~=1)
  error('Feature type is not a_cnch_cnch3_va_vna_va3_vna3');
end

mito_dir = [get_reconstruction_dir(config), mitochondria_config.dir];

if(mitochondria_config.is_verbose)
  fprintf('\nSTART: mitochondria_segment_apply_classifier\n');
end

if(mitochondria_config.is_verbose)
  fprintf('Loading segment features ... ');
end
features_stack=[];
for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  if(mitochondria_config.is_verbose)
    fprintf('plane:%d\n', case_id);
  end
  image_prefixes = get_image_prefixes_subdirs(config, case_id);
  image_prefix = image_prefixes{1};

  load_file_name = [mito_dir, image_prefix, ...
    '.mito.seg', segment_config.mito_seg_suffix, ...
    '.features', feature_config.features_suffix, ...
    '.mat'];
  segment_features = load2(load_file_name, 'features');
  features_stack=[features_stack; i*ones(length(segment_features.features),1), ...
    segment_features.features ];

end
if(mitochondria_config.is_verbose)
  fprintf('done.\n');
end

if(mitochondria_config.is_verbose)
  fprintf('Loading model  ... ');
end

load_file_name = [mito_dir, ...
  '.mito.seg', segment_config.mito_seg_suffix, ...
  '.features', feature_config.features_suffix, ...
  '.gt', segment_config.gt_suffix, ...
  '.model', model_config.model_suffix, ...
  '.mat'];
model = load2(load_file_name, 'model_config');

if(mitochondria_config.is_verbose)
  fprintf('done.\n');
end

if(mitochondria_config.is_verbose)
  fprintf('Apply classifier ... ');
end
samples_test=[features_stack(:,3:end)];
boosting_stack=Classify(model.model_config.GLearners, model.model_config.GWeights, samples_test')';

if(mitochondria_config.is_verbose)
  fprintf('done.\n');
end

if(mitochondria_config.is_verbose)
  fprintf('Saving confidence values ... ');
end

for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  image_prefixes = get_image_prefixes_subdirs(config, case_id);
  image_prefix = image_prefixes{1};
  seg_mt_conf_vals=boosting_stack(find(features_stack(:,1)==i));

  save_file_name = [mito_dir, image_prefix,...
    '.mito.seg', segment_config.mito_seg_suffix, ...
    '.features', feature_config.features_suffix, ...
    '.gt', segment_config.gt_suffix, ...
    '.model', model_config.model_suffix, ...
    '.conf_values',...
    '.mat'];
  save2(save_file_name, 'seg_mt_conf_vals');

end

if(mitochondria_config.is_verbose)
  fprintf('done.\n');
end

if(mitochondria_config.is_verbose)
  fprintf('\nSTOP: mitochondria_segment_apply_classifier\n');
end
return;
end
