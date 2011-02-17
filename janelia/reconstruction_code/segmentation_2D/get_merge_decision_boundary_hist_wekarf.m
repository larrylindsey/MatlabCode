function [label_map_new, sp_to_seg_map] = ...
  get_merge_decision_boundary_hist_wekarf(label_map, boundary_map, ...
  boundary_sets, merge_criterion_param, image, config) %#ok<INUSD,INUSL>
% to_merge_sets = get_merge_decision_boundary_hist_boost(label_map, boundary_map, ...
%   boundary_sets, merge_criterion_param)
%
global config_global

% seg_model_weka = load('/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/scripts/linkage_experiments/classification_exp/seg_wekarf_boundary_hist_300_1_ex0.3_596_605.mat');
seg_model_weka = load('/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/scripts/linkage_experiments/classification_exp/seg_wekarf_boundary_hist_300_1_ex0.3_l30_596_605.mat');

interior_map = {};
interior_map{1} = 1 - boundary_map;
% interior_map{end+1} = imrotate(filter_image2(imrotate(image, 90), ...
%   'db_nwo700_cbf2_LDA7_mf7_ps1\_0.25_abf_neg', config), -90);
% interior_map{end+1} = imrotate(filter_image2(imrotate(image, 180), ...
%   'db_nwo700_cbf2_LDA7_mf7_ps1\_0.25_abf_neg', config), -180);
% interior_map{end+1} = imrotate(filter_image2(imrotate(image, -90), ...
%   'db_nwo700_cbf2_LDA7_mf7_ps1\_0.25_abf_neg', config), 90);
% interior_map{end+1} = ones(size(boundary_map));
% interior_map{end}(2:end, 2:end) = 1 - boundary_map(1:end-1, 1:end-1);
% interior_map{end+1} = ones(size(boundary_map));
% interior_map{end}(1:end-1, 1:end-1) = 1 - boundary_map(2:end, 2:end);
% interior_map{end+1} = ones(size(boundary_map));
% interior_map{end}(1:end-1, 2:end) = 1 - boundary_map(2:end, 1:end-1);
% interior_map{end+1} = ones(size(boundary_map));
% interior_map{end}(2:end, 1:end-1) = 1 - boundary_map(1:end-1, 2:end);
% interior_map{end+1} = ones(size(boundary_map));
% interior_map{end}(3:end, 3:end) = 1 - boundary_map(1:end-2, 1:end-2);
% interior_map{end+1} = ones(size(boundary_map));
% interior_map{end}(1:end-2, 1:end-2) = 1 - boundary_map(3:end, 3:end);
% interior_map{end+1} = ones(size(boundary_map));
% interior_map{end}(1:end-2, 3:end) = 1 - boundary_map(3:end, 1:end-2);
% interior_map{end+1} = ones(size(boundary_map));
% interior_map{end}(3:end, 1:end-2) = 1 - boundary_map(1:end-2, 3:end);
% interior_map{end+1} = 1 - boundary_map;
% interior_map{end}(interior_map{end}<=0.99) = ...
%   interior_map{end}(interior_map{end}<=0.99) + 0.01;
% interior_map{end+1} = 1 - boundary_map;
% interior_map{end}(interior_map{end}>=0.01) = ...
%   interior_map{end}(interior_map{end}>=0.01) - 0.01;
% interior_map{end+1} = 1 - boundary_map + 0.1*(rand(size(boundary_map))-0.5);
% interior_map{end}(interior_map{end}>1) = 1;
% interior_map{end}(interior_map{end}<0) = 0;

n_iteration = 0;
no_change = false;
label_map_new = label_map;
sp_ids = [0; nonzeros(unique(label_map(:)))];
sp_to_seg_map = [sp_ids, sp_ids];
while(~no_change)
  n_iteration = n_iteration + 1;
  to_merge_sets = [];
  for i = 1:length(interior_map)
    boundary_hist = collect_segment_pair_stats_boundary_hist(...
      uint32(label_map_new), uint8(255*interior_map{i}), ...
      uint8(10));
    boundary_hist = boundary_hist(min(boundary_hist(:,1:2), [], 2)>0, :);
    boundary_hist = boundary_hist(boundary_hist(:,1)~=boundary_hist(:,2), :);
    features_boundary = boundary_hist(:, 3:end);
    
    features = features_boundary;
%     if(size(features,2)~=model.n_dimension_insection)
%       error('Segment merge feature dimensions don''t match during training and testing.');
%     end
    write_arff('seg_temp', features, ones(size(features,1),1), ...
      [config_global.temp_dir, 'seg_temp.arff']);
    data_weka = loadARFF([config_global.temp_dir, 'seg_temp.arff']);
    [junk, confidence_boundary] = wekaClassify(data_weka, seg_model_weka.classifier); %#ok<ASGLU>
    
    mean_boundary_test = sum(features .* repmat((0:10:255)/255, ...
      [size(features,1), 1]), 2)./sum(features, 2);
    to_force = mean_boundary_test<0.3 | sum(features,2)<=30;
    confidence_boundary(to_force, 1) = 0;
    confidence_boundary(to_force, 2) = 1;
    
    to_merge_sets = [to_merge_sets; boundary_hist(...
      confidence_boundary(:,1)>merge_criterion_param.threshold, [1 2])]; %#ok<AGROW>
  end
  
  if(~isempty(to_merge_sets))
    % make mergers
    max_label = max(label_map_new(:));
    % merger adjacency matrix for connected components analysis
    A = sparse([], [], [], max_label, max_label, 2*size(to_merge_sets,1));
    A(sub2ind(size(A), to_merge_sets(:,1), to_merge_sets(:,2))) = 1; %#ok<SPRIX>
    A(sub2ind(size(A), to_merge_sets(:,2), to_merge_sets(:,1))) = 1; %#ok<SPRIX>
    
    seg_label_mapping = [0; get_connected_components_c(A, 0.5)];
    
    % generate new label map
    label_map_new = seg_label_mapping(1+label_map_new);
    label_map_new = double(remove_merged_boundaries_2D(uint32(label_map_new)));
    
    % generate superpixel to segment map
    sp_to_seg_map(:,2) = seg_label_mapping(1+sp_to_seg_map(:,2));
  else
    no_change = true;
  end
  if(n_iteration>20)
    break;
  end
end

return
end
