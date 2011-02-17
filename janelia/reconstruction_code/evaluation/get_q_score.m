function get_q_score(config)
% get_q_score(config)
% Evaluate reconstructed volume using q_score
%
% Takes the reconstruction created in seg and evaluates using q_score.
% Calls qscore3.m - Yuriy Mishchenko
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04142008  init code
%

fprintf('Evaluating the reconstructed volume with Q Score ..\n');

linkage_config = config.linkage;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
apply_config = linkage_config.apply;
volume_config = config.final_volume;
evaluate_config = config.evaluate_volume;

load2([get_stack_dir(config), evaluate_config.annotations_dir, evaluate_config.groundtruth_file], 'seg');
seg_gt = seg;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(volume_config.dir, 'dir')~=7)
  mkdir2(volume_config.dir);
end;
cd(prev_dir);

t=1;
fprintf('%d ', t);
load2([get_reconstruction_dir(config), volume_config.dir, 'seg', '.',  model_config.type, '.', feature_config.type, ...
  apply_config.model_suffix, '.lk', num2str(volume_config.link_thresholds(t)), ...
  apply_config.segmentation_suffix, '.mat'], 'seg');

[q_score, q_score_d] = qscore3(seg, seg_gt, false);
for t = 2:length(volume_config.link_thresholds)
  fprintf('%d ', t);
  if(mod(t,20)==0)
    fprintf('\n');
  end;
  
  load2([get_reconstruction_dir(config), volume_config.dir, '/', 'seg', '.',  model_config.type, '.', feature_config.type, ...
    apply_config.model_suffix, '.lk', num2str(volume_config.link_thresholds(t)), ...
    apply_config.segmentation_suffix, '.mat'], 'seg');

  [q, d] = qscore3(seg, seg_gt, false);
  q_score(end+1) = q;
  q_score_d(end+1) = d;

end;

link_thresholds = volume_config.link_thresholds;
save2([get_reconstruction_dir(config), volume_config.dir, 'q_score', '.',  model_config.type, '.', feature_config.type, ...
  apply_config.model_suffix, apply_config.segmentation_suffix, '.mat'], 'q_score', 'q_score_d', 'link_thresholds');

fprintf('\ndone.\n');


return
end