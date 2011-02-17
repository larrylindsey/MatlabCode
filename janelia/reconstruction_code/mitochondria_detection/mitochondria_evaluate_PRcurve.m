function mitochondria_evaluate_PRcurve(config)
% mitochondria_evaluate_PRcurve(config)
% Generate PR curve
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


error('Not yet (March 17th 2009) implemented.');

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
  fprintf('\nSTART: mitochondria_evaluate_PRcurve\n');
end

avg.precision = zeros(length(segment_config.PR_curve_suffixes),1);
avg.recall = zeros(length(segment_config.PR_curve_suffixes),1);

for suffix_num = 1:length(segment_config.PR_curve_suffixes);
    if(mitochondria_config.is_verbose)
        fprintf('For parameters:%d\n', suffix_num);
    end
  
precision=zeros(length(stack_config.case_ids),1);
recall=zeros(length(stack_config.case_ids),1);

for i = 1:length(stack_config.case_ids)
    case_id = stack_config.case_ids(i);
    if(mitochondria_config.is_verbose)
        fprintf('plane:%d\n', case_id);
    end
    image_prefixes = get_image_prefixes_subdirs(config, case_id);
    image_prefix = image_prefixes{1};

load_file_name = [mito_dir, image_prefix,...
  segment_config.PR_curve_suffixes, ...
  '.mat'];
variable_name = horzcat('precision',segment_config.PR_curve_error);
x = load2(load_file_name, variable_name );
variable_name = horzcat('recall',segment_config.PR_curve_error);
y = load2(load_file_name, variable_name);

  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end
  
end

avg.precision(suffix_num) = mean(precision);
avg.recall(suffix_num) = mean(recall);
end

avg.precision
avg.recall

  if(mitochondria_config.is_verbose)
  fprintf('\nSTOP: mitochondria_evaluate_PRcurve\n');
end
return;
end
  
  