function mitochondria_segment_label_MRF(config)
% mitochondria_segment_label_MRF(config)
% Uses Markov Random Field to classify segments based on linkage weights
% and mitochondria confidence values
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
  fprintf('\nSTART: mitochondria_segment_label_MRF\n');
end

if(mitochondria_config.is_verbose)
  fprintf('Loading segment confidence values and weights ... ');
end
cur_tot_number_segs=0;
prev_tot_number_segs=0;
seg_conf_values=[];
seg_linkage_weights=[];

for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  if(mitochondria_config.is_verbose)
    fprintf('plane:%d\n', case_id);
  end
  image_prefixes = get_image_prefixes_subdirs(config, case_id);
  image_prefix = image_prefixes{1};

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
  seg_conf_values=[seg_conf_values;seg_conf.seg_mt_conf_vals];

  %%%% load segment linkage weights for image
  load_file_name = [mito_dir, image_prefix, ...
    '.mito.seg', segment_config.mito_seg_suffix, ...
    '.linkage_conf', ...
    '.mat'];
  seg_links = load2(load_file_name, 'linkage_features');

  seg_links.linkage_features(:,1)=seg_links.linkage_features(:,1)+cur_tot_number_segs;
  seg_links.linkage_features(:,2)=seg_links.linkage_features(:,2)+(cur_tot_number_segs+num_segs)*seg_links.linkage_features(:,3);
  seg_links.linkage_features(:,2)=seg_links.linkage_features(:,2)+(prev_tot_number_segs)*(1-seg_links.linkage_features(:,3));
  seg_links.linkage_features(:,3)=[];
  seg_linkage_weights=[seg_linkage_weights;seg_links.linkage_features];

  prev_tot_number_segs=cur_tot_number_segs;
  cur_tot_number_segs=cur_tot_number_segs+num_segs;
end
if(mitochondria_config.is_verbose)
  fprintf('done.\n');
end


if(mitochondria_config.is_verbose)
  fprintf('Applying Markov Random Field ... ');
end

%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
%seg_conf_thresh=-1;
%avg_linkage_weight=150;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seg_conf_values=seg_conf_values-mrf_config.seg_conf_thresh;
seg_linkage_weights=seg_linkage_weights(seg_linkage_weights(:,5)>=.1,:);

%%
label_prior = repmat(abs(seg_conf_values), [1 2]);
label_prior(seg_conf_values<0,1) = 0;
label_prior(seg_conf_values>=0,2) = 0;
label_prior = int32(label_prior*100);

%%
pairwise_smooth = seg_linkage_weights(:, [1 2 5]);
pairwise_smooth(:,1:2) = pairwise_smooth(:,1:2) - 1;
pairwise_smooth(:,3) = mrf_config.avg_linkage_weight*pairwise_smooth(:,3);
pairwise_smooth = int32(pairwise_smooth);

label = graphcut_mex(label_prior, pairwise_smooth, 1);

if(mitochondria_config.is_verbose)
  fprintf('done.\n');
end

if(mitochondria_config.is_verbose)
  fprintf('Saving mitochondria segments ... ');
end
cur_tot_number_segs=0;
%%% Find segments labeled as mitochondria
for i = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  if(mitochondria_config.is_verbose)
    fprintf('plane:%d\n', case_id);
  end
  image_prefixes = get_image_prefixes_subdirs(config, case_id);
  image_prefix = image_prefixes{1};

  load_file_name = [mito_dir, image_prefix, ...
    '.mito.seg', segment_config.mito_seg_suffix,...
    '.mat'];
  segments = load2(load_file_name, 'label_map');
  num_segs=max(segments.label_map(:));
  indices=[cur_tot_number_segs+1:cur_tot_number_segs+num_segs];
  seg_labels=label(indices);
  cur_tot_number_segs=cur_tot_number_segs+num_segs;

  %% Save mt segments
  save_file_name = [mito_dir, image_prefix,...
    '.mito.seg', segment_config.mito_seg_suffix, ...
    '.features', feature_config.features_suffix, ...
    '.gt', segment_config.gt_suffix, ...
    '.model', model_config.model_suffix, ...
    '.mrf.', num2str(mrf_config.seg_conf_thresh), ...
    '.', num2str(mrf_config.avg_linkage_weight), ...
    '.mat'];
  save2(save_file_name, 'seg_labels');
end
if(mitochondria_config.is_verbose)
  fprintf('done.\n');
end

if(mitochondria_config.is_verbose)
  fprintf('\nSTOP: mitochondria_segment_label_MRF\n');
end
return;
end
