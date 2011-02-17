function get_rand_score(config)
% get_rand_score(config)
% Evaluate reconstructed volume using Rand
%
% Takes the reconstruction created in seg and evaluates using Rand.
% Calls partitionDistance.m - Alex Genkin
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04142008  init code
%

fprintf('Evaluating the reconstructed volume with Rand Score ..\n');

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

[a_not_b, b_not_a] = partitionDist(seg, seg_gt);
[wRand_oversegm,wRand_undersegm] = wRand(seg, seg_gt);
for t = 2:length(volume_config.link_thresholds)
  fprintf('%d ', t);
  if(mod(t,20)==0)
    fprintf('\n');
  end;
  
  load2([get_reconstruction_dir(config), volume_config.dir, '/', 'seg', '.',  model_config.type, '.', feature_config.type, ...
    apply_config.model_suffix, '.lk', num2str(volume_config.link_thresholds(t)), ...
    apply_config.segmentation_suffix, '.mat'], 'seg');

  [a_n_b, b_n_a] = partitionDist(seg, seg_gt);
  a_not_b(end+1) = a_n_b;
  b_not_a(end+1) = b_n_a;

  %weighted Rand
  [overs,unders] = wRand(seg, seg_gt);
  wRand_oversegm(end+1) = overs;
  wRand_undersegm(end+1) = unders;

end;

link_thresholds = volume_config.link_thresholds;
save2([get_reconstruction_dir(config), volume_config.dir, 'rand_score', '.',  model_config.type, '.', feature_config.type, ...
  apply_config.model_suffix, apply_config.segmentation_suffix, '.mat'],...
  'a_not_b', 'b_not_a', 'wRand_oversegm','wRand_undersegm', 'link_thresholds');

% display
[r_best,r_best_idx] = min(a_not_b+b_not_a);
fprintf('\n\n');disp(['Best Rand ' num2str(r_best) ' at threshold ' num2str(link_thresholds(r_best_idx))]);
[wr_best,wr_best_idx] = min( wRand_oversegm + wRand_undersegm );
fprintf('\n\n');disp(['Best Weighted Rand ' num2str(wr_best) ' at threshold ' num2str(link_thresholds(wr_best_idx))]);

figure('name',['Rand: ' config.stack.name]),
subplot(1,2,1),plot(link_thresholds',[a_not_b',b_not_a',a_not_b'+b_not_a'])
title('Rand'),legend('undersegm','oversegm','Rand')
subplot(1,2,2),plot(link_thresholds', [wRand_undersegm', wRand_oversegm', wRand_oversegm'+wRand_undersegm'])
title('Weighted Rand'),legend('undersegm','oversegm','wRand')

fprintf('\ndone.\n');


return
end