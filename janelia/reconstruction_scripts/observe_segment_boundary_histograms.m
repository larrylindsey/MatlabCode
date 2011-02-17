load ../lamina.OTO.3500x.zhiyuan.May232008.3_5/train_segmentation/crop.a.003seg_boundary_sets_gt.v0.mat
% load ../lamina.OTO.3500x.zhiyuan.May232008.3_5/train_segmentation/crop.a.004seg_boundary_sets_gt.v0.mat
% load ../lamina.OTO.3500x.zhiyuan.May232008.3_5/train_segmentation/crop.a.005seg_boundary_sets_gt.v0.mat

% load ../lamina.OTO.3500x.zhiyuan.May232008.3_5/train_segmentation/crop.a.003.seg_boundary_sets_gt.v_heq_e2.mat
% load ../lamina.OTO.3500x.zhiyuan.May232008.3_5/train_segmentation/crop.a.004.seg_boundary_sets_gt.v_heq_e2.mat
% load ../lamina.OTO.3500x.zhiyuan.May232008.3_5/train_segmentation/crop.a.005.seg_boundary_sets_gt.v_heq_e2.mat

figure(13);
intensity_bins = 0:0.1:1;
h = hist(seg_boundary_sets_gt(1).boundary_set, intensity_bins);
plot(intensity_bins, h/sum(h));
for i = 2:length(seg_boundary_sets_gt)
  if(seg_boundary_sets_gt(i).segment_label_0==1045 || ...
      seg_boundary_sets_gt(i).segment_label_1==1045)
    continue;
  end
  if(size(seg_boundary_sets_gt(i).boundary_set, 2)<20)
    continue;
  end
  h = hist(seg_boundary_sets_gt(i).boundary_set, intensity_bins);
  hold on;
  plot(intensity_bins, h/sum(h));
  hold off;
end

%%
train_seg_dir = '../lamina.OTO.3500x.zhiyuan.May232008.3_5/train_segmentation/';
% load([train_seg_dir, ...
%   'crop.a.003.seg_boundary_sets.v_heq_e2.gs_l_sp_T0.5_L150.gs_l_T0.2_L100.v_heq_e1_m0_d5_a20.v_heq_e2.mat']);
% load([train_seg_dir, ...
%   'crop.a.003.seg_boundary_sets.v_heq.gs_l_sp_T0.5_L150.gs_l_T0.2_L100.v_heq_e1_m0_d5_a20.v_heq_e2.mat']);
load([train_seg_dir, ...
  'crop.a.003.seg_boundary_sets.v0.gs_l_sp_T0.5_L150.gs_l_T0.2_L100.v_heq_e1_m0_d5_a20.v_heq_e2.mat']);

figure;
intensity_bins = 0:0.1:1;
h = hist(seg_boundary_sets.matched(1).boundary_set, intensity_bins);
plot(intensity_bins, h/sum(h));
boundary_means_matched = [];
for i = 2:length(seg_boundary_sets.matched)
  if(seg_boundary_sets.matched(i).segment_label_0==1045 || ...
      seg_boundary_sets.matched(i).segment_label_1==1045)
    continue;
  end
  if(size(seg_boundary_sets.matched(i).boundary_set, 2)<20)
    continue;
  end
  h = hist(seg_boundary_sets.matched(i).boundary_set(3,:), intensity_bins);
  hold on;
  plot(intensity_bins, h/sum(h));
  hold off;
  
  boundary_means_matched(end+1) = mean(seg_boundary_sets.matched(i).boundary_set(3,:));
end
figure; plot(intensity_bins, hist(boundary_means_matched, intensity_bins));

figure;
intensity_bins = 0:0.1:1;
h = hist(seg_boundary_sets.unmatched(1).boundary_set, intensity_bins);
plot(intensity_bins, h/sum(h));
boundary_means_unmatched = [];
for i = 2:length(seg_boundary_sets.unmatched)
  if(seg_boundary_sets.unmatched(i).segment_label_0==1045 || ...
      seg_boundary_sets.unmatched(i).segment_label_1==1045)
    continue;
  end
  if(size(seg_boundary_sets.unmatched(i).boundary_set, 2)<20)
    continue;
  end
  h = hist(seg_boundary_sets.unmatched(i).boundary_set(3,:), intensity_bins);
  hold on;
  plot(intensity_bins, h/sum(h)); 
  hold off;
  
  boundary_means_unmatched(end+1) = mean(seg_boundary_sets.unmatched(i).boundary_set(3,:));  
end
figure; plot(intensity_bins, hist(boundary_means_unmatched, intensity_bins));

%%
train_seg_dir = '../lamina.OTO.3500x.zhiyuan.May232008.3_5/train_segmentation/';
seg_id = {'0.1_0_b0','0.21_0.2_b0','0.31_0.3_b0','0.41_0.4_b0','0.51_0.5_b0','0.61_0.6_b0', '0.71_0.7_b25'};
min_boundary_lengths = [6, 6, 6, 6, 6, 6, 25];
plot_legend = {'b.', 'r.', 'g.', 'm.', 'c.', 'y.', 'k.'};
for t = 1:length(seg_id)
%   load([train_seg_dir, ...
%   'crop.a.003.seg_boundary_sets.v0.gs_l_sp_T0.1_L150.gs_amb_T', seg_id{t}, ...
%   '.v_heq_o1_m0_d5_a40.v_heq.v_heq_e2.mat']);
  load([train_seg_dir, ...
  'crop.a.003.seg_boundary_sets.v0.gs_amb_T', seg_id{t}, ...
  '.v_heq_o1_m0_d5_a40.v_heq.mat']);

  figure(100*t);
  intensity_bins = 0:0.1:1;
  h = hist(seg_boundary_sets.matched(1).boundary_set, intensity_bins);
  plot(intensity_bins, cumsum(h/sum(h)));
  boundary_means_matched = [];
  boundary_std_matched = [];
  for i = 2:length(seg_boundary_sets.matched)
    if(seg_boundary_sets.matched(i).segment_label_0==1045 || ...
        seg_boundary_sets.matched(i).segment_label_1==1045)
      continue;
    end
    if(seg_boundary_sets.matched(i).segment_label_0==...
        seg_boundary_sets.matched(i).segment_label_1)
      continue;
    end
    if(size(seg_boundary_sets.matched(i).boundary_set, 2)<min_boundary_lengths(t))
      continue;
    end
    h = hist(seg_boundary_sets.matched(i).boundary_set(3,:), intensity_bins);
    hold on;
    plot(intensity_bins, cumsum(h/sum(h)));
    hold off;

    boundary_means_matched(end+1,:) = [size(seg_boundary_sets.matched(i).boundary_set, 2), ...
      mean(seg_boundary_sets.matched(i).boundary_set(3,:))];
    boundary_std_matched(end+1,:) = [size(seg_boundary_sets.matched(i).boundary_set, 2), ...
      std(seg_boundary_sets.matched(i).boundary_set(3,:))];
  end
  title(['Matched ', seg_id{t}]);
  axis([0 1 0 1]);
  h = hist(boundary_means_matched(:,2), intensity_bins);
  figure(1); plot(intensity_bins, cumsum(h/sum(h)));
  hold on;
  title(['Matched means ', seg_id{t}]);
  axis([0 1 0 1]);
  h = hist(boundary_std_matched(:,2), intensity_bins);
  figure(2); plot(intensity_bins, cumsum(h/sum(h)));
  hold on;
  title(['Matched stds ', seg_id{t}]);
  axis([0 1 0 1]);
  figure(11); plot(boundary_means_matched(:,1), boundary_means_matched(:,2), ...
    plot_legend{t});
  hold on;
  title(['Matched means vs. boundary lengths', seg_id{t}]);
  figure(12); plot(boundary_std_matched(:,1), boundary_std_matched(:,2), ...
    plot_legend{t});
  hold on;
  title(['Matched stds vs. boundary lengths', seg_id{t}]);

  figure(100*t+1);
  intensity_bins = 0:0.1:1;
  h = hist(seg_boundary_sets.unmatched(1).boundary_set, intensity_bins);
  plot(intensity_bins, cumsum(h/sum(h)));
  boundary_means_unmatched = [];
  boundary_std_unmatched = [];
  for i = 2:length(seg_boundary_sets.unmatched)
    if(seg_boundary_sets.unmatched(i).segment_label_0==1045 || ...
        seg_boundary_sets.unmatched(i).segment_label_1==1045)
      continue;
    end
    if(seg_boundary_sets.unmatched(i).segment_label_0==...
        seg_boundary_sets.unmatched(i).segment_label_1)
      continue;
    end
    if(size(seg_boundary_sets.unmatched(i).boundary_set, 2)<min_boundary_lengths(t))
      continue;
    end
    h = hist(seg_boundary_sets.unmatched(i).boundary_set(3,:), intensity_bins);
    hold on;
    plot(intensity_bins, cumsum(h/sum(h)));
    hold off;

    boundary_means_unmatched(end+1,:) = [size(seg_boundary_sets.unmatched(i).boundary_set, 2), ...
      mean(seg_boundary_sets.unmatched(i).boundary_set(3,:))];
    boundary_std_unmatched(end+1,:) = [size(seg_boundary_sets.unmatched(i).boundary_set, 2), ...
    std(seg_boundary_sets.unmatched(i).boundary_set(3,:))];
  end
  title(['Unmatched ', seg_id{t}]);
  axis([0 1 0 1]);
  h = hist(boundary_means_unmatched(:,2), intensity_bins);
  figure(3); plot(intensity_bins, cumsum(h/sum(h)));
  hold on
  title(['Unmatched means ', seg_id{t}]);
  axis([0 1 0 1]);
  h = hist(boundary_std_unmatched(:,2), intensity_bins);
  figure(4); plot(intensity_bins, cumsum(h/sum(h)));
  hold on
  title(['Unmatched stds ', seg_id{t}]);
  axis([0 1 0 1]);
  figure(13); plot(boundary_means_unmatched(:,1), boundary_means_unmatched(:,2), ...
    plot_legend{t});
  hold on;
  title(['Unmatched means vs. boundary lengths', seg_id{t}]);
  figure(14); plot(boundary_std_unmatched(:,1), boundary_std_unmatched(:,2), ...
    plot_legend{t});
  title(['Unmatched stds vs. boundary lengths', seg_id{t}]);
  hold on;
end

figure(1); hold off;
figure(2); hold off;
figure(3); hold off;
figure(4); hold off;
figure(11); hold off;
figure(12); hold off;
figure(13); hold off;
figure(14); hold off;

%%
train_seg_dir = '../lamina.OTO.3500x.zhiyuan.May232008.3_5/train_segmentation/';
seg_id = {'0.21_0.2_b0','0.31_0.3_b0','0.41_0.4_b0','0.51_0.5_b0','0.61_0.6_b0', '0.71_0.7_b25'};

length_bins = 1:5:100;
for t = 1:length(seg_id)
%   load([train_seg_dir, ...
%   'crop.a.003.seg_boundary_sets.v_heq.gs_l_sp_T0.1_L150.gs_amb_T', seg_id{t}, ...
%   '.v_heq_o1_m0_d5_a40.v_heq.v_heq_e2.mat']);
  load([train_seg_dir, ...
  'crop.a.003.seg_boundary_sets.v0.gs_amb_T', seg_id{t}, ...
  '.v_heq_o1_m0_d5_a40.v_heq.mat']);

  boundary_lengths_matched = [];
  for i = 2:length(seg_boundary_sets.matched)
    if(seg_boundary_sets.matched(i).segment_label_0==1045 || ...
        seg_boundary_sets.matched(i).segment_label_1==1045)
      continue;
    end
    if(seg_boundary_sets.matched(i).segment_label_0==...
        seg_boundary_sets.matched(i).segment_label_1)
      continue;
    end
    size(seg_boundary_sets.matched(i).boundary_set, 2)
    boundary_lengths_matched(end+1) = size(seg_boundary_sets.matched(i).boundary_set, 2);
  end
  h = hist(boundary_lengths_matched, length_bins);
  figure(11); plot(length_bins, h/sum(h));
  hold on;
  title(['Matched boundary lengths ', seg_id{t}]);

  boundary_lengths_unmatched = [];
  for i = 2:length(seg_boundary_sets.unmatched)
    if(seg_boundary_sets.unmatched(i).segment_label_0==1045 || ...
        seg_boundary_sets.unmatched(i).segment_label_1==1045)
      continue;
    end
    if(seg_boundary_sets.unmatched(i).segment_label_0==...
        seg_boundary_sets.unmatched(i).segment_label_1)
      continue;
    end
    boundary_lengths_unmatched(end+1) = size(seg_boundary_sets.unmatched(i).boundary_set, 2);
  end
  h = hist(boundary_lengths_unmatched, length_bins);
  figure(12); plot(length_bins, h/sum(h));
  hold on
  title(['Unmatched boundary lengths ', seg_id{t}]);
end

figure(11); hold off;
figure(12); hold off;
