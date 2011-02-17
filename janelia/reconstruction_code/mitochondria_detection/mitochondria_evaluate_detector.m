function mitochondria_evaluate_detector(config)
% mitochondria_evaluate_detector(config)
% Evaluate mitochondria detection against manual annotation
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  05052008  init code
%

stack_config = config.stack;
mitochondria_config = config.mitochondria;
train_config = mitochondria_config.train;
evaluate_config = mitochondria_config.evaluate;

image_dir = [get_reconstruction_dir(config), mitochondria_config.dir, train_config.dir];

if(isfield(evaluate_config, 'verbose'))
  is_verbose = evaluate_config.verbose;
else
  is_verbose = true;
end;

THRESHOLDS = -10:0.1:10; % ;
n_threshold = length(THRESHOLDS);

precision = zeros(1, n_threshold);
total_prec = zeros(1, n_threshold);
recall = zeros(1, n_threshold);
total_recl = zeros(1, n_threshold);

for case_id = stack_config.case_ids
  fprintf('case %d: ', case_id);
  
  load2(sprintf([get_reconstruction_dir(config), mitochondria_config.dir, stack_config.image_prefix, ...
    evaluate_config.suffix, '.mat'], case_id), 'mitochondria_detection_confidence');
  
  load2( sprintf([image_dir, train_config.image_prefix, train_config.annot_suffix, '.mat'], ...
    case_id), 'zone_mask');
  zone_mask = zone_mask>0;

  for t = 1:n_threshold
    m = mitochondria_detection_confidence>THRESHOLDS(t);
    m = bwareaopen(m, round(625));
    
    precision(t) = precision(t) + sum(m(:) & zone_mask(:));
    total_prec(t) = total_prec(t) + sum(m(:));
    
    recall(t) = recall(t) + sum(m(:) & zone_mask(:));
    total_recl(t) = total_recl(t) + sum(zone_mask(:));
  end

%   m_low_conf = mitochondria_detection_confidence>0;
%   label_low_conf = bwlabel(m_low_conf);
%   
%   m_high_conf = mitochondria_detection_confidence>5;
%   label_high_conf = bwlabel(m_high_conf);
%   stats_high_conf = regionprops(label_high_conf, 'Area');
%   id_high_conf  = find([stats_high_conf(:).Area]>1000);
%   seeds_high_conf = ismember(label_high_conf, id_high_conf);
%   
%   seeded_label_id = unique(label_low_conf(seeds_high_conf));
%   
%   m_filtered = ismember(label_low_conf, seeded_label_id);
%   
%   precision = precision + sum(m_filtered(:) & zone_mask(:));
%   total_prec = total_prec + sum(m_filtered(:));
% 
%   recall = recall + sum(m_filtered(:) & zone_mask(:));
%   total_recl = total_recl + sum(zone_mask(:));
  
  fprintf('\n');
end

if(is_verbose)
  P = precision./total_prec;
  R = recall./total_recl;

  figure(500); hold on;
  plot(R, P, 'ro-');
  axis([0 1 0 1]);
  hold off;
  title('Mitochondria detection:  Precision-Recall curve');
end
  
return;
end
