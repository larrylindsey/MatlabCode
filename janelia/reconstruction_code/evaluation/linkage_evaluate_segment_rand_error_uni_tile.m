function linkage_evaluate_segment_rand_error_uni_tile(config)
% linkage_evaluate_segment_rand_error(config)
% Train and test on the training data using k-fold protocol - (k-1)/k
% fraction of training data's sections are used for training and rest for
% testing.
%
% Pratim Ghosh,
% Summer intern, Janelia Farm Research Campus, HHMI.
% Univ. of California, Santa Barbara.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v1  07092008  init code
% v2  06222009  split data section-wise
% v3  07282009  include modified segmentations.
%

config_0 = config;
stack_config = config.stack;
linkage_config = config.linkage;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
train_config = linkage_config.train;
evaluate_config = linkage_config.evaluate;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(linkage_config.dir, 'dir')~=7)
  mkdir2(linkage_config.dir);
end;
cd(prev_dir);

fprintf('Loading manual annotation .. ');
manual_annotation = load2(train_config.manual_annotation_file);
fprintf('done\n');

if(isfield(evaluate_config, 'segmentation_file') && ...
    ~isempty(evaluate_config.segmentation_file))
  fprintf('Loading segmentation to be used for training and testing .. ');
  segmentation = load2(evaluate_config.segmentation_file);
  fprintf('done\n');
else
  segmentation = manual_annotation;
end

if(isfield(evaluate_config, 'original_segmentation_file') && ...
    ~isempty(evaluate_config.original_segmentation_file))
  fprintf('Loading original segmentation to be used for assignment to gt ... ');
  original_segmentation = load2(evaluate_config.original_segmentation_file);
  fprintf('done\n');
else
  original_segmentation = segmentation;
end

n_section = length(manual_annotation.al);

% iteration over the k folds
config.stack.roi = [];
config.stack.align.is_prealigned = true;
config.stack.dir = evaluate_config.temporary_dir;

% save images in a temporary directory
prev_dir = pwd2;
cd(config.stack.dir);
if(exist(stack_config.name, 'dir')~=7)
  mkdir2(stack_config.name);
end;
cd(prev_dir);
config.stack.case_ids = 1:n_section;
config.segmentation_2D.method = 'link_eval_k';
config.superpixel_2_seg.method = 'link_eval_k';
config.superpixel_choose.choice.method = 'link_eval_k';
config.superpixel_choose.choice.seg_suffix = '.eval_k';
config.segmentation_choose.choice.method = 'link_eval_k';
config.segmentation_choose.choice.seg_suffix = '.eval_k';
config.linkage.train.save_suffix = [config_0.linkage.train.save_suffix, 'k'];
config.linkage.apply.model_suffix = config.linkage.train.save_suffix;
% config.linkage.apply.segmentation_suffix = evaluate_config.segmentation_suffix;

image_dir = [config.stack.dir, stack_config.name, '/'];
seg_dir = get_segmentation_dir(config);
fprintf('Assigning automatic segmentation to ground-truth segments ...\n');
label_mapping_as_to_gt = {};
for z = 1:n_section
  fprintf('z: %d\n', z);
  imwrite(manual_annotation.al{z}, sprintf([image_dir, stack_config.image_prefix, ...
    stack_config.image_suffix], z));
  
  % Assign automatic segmentation map's segments to ground-truth segments.
  % If an A.S. segment overlaps substantially with more than one
  % ground-truth segment then disregard it during training and testing.
  [label_pairs, overlap_area] = count_row_occurence(...
    [original_segmentation.wcat{z}(:), manual_annotation.wcat{z}(:)]);
  overlap_area = overlap_area(label_pairs(:,1)>0 & label_pairs(:,2)>0);
  label_pairs = label_pairs(label_pairs(:,1)>0 & label_pairs(:,2)>0, :);
  
  [junk, a_1] = count_row_occurence(original_segmentation.wcat{z}(:));
  area_1 = [];
  area_1(1+junk) = a_1; %#ok<AGROW>
  [junk, a_2] = count_row_occurence(manual_annotation.wcat{z}(:));
  area_2 = [];
  area_2(1+junk) = a_2; %#ok<AGROW>
  mapping = correspond_labels_many_to_one(label_pairs(:,1), ...
    label_pairs(:,2), overlap_area, area_1', area_2', 200, 0.1);
  label_mapping_as_to_gt{z} = zeros(max(original_segmentation.wcat{z}(:))+1, 1); %#ok<AGROW>
  label_mapping_as_to_gt{z}(mapping(:,1) + 1) = mapping(:,2); %#ok<AGROW>
  
  is_mapped_to_0 = find(label_mapping_as_to_gt{z}==0) - 1;
  segmentation.wcat{z}(ismember(segmentation.wcat{z}, is_mapped_to_0)) = 0;
  
  label_map = segmentation.wcat{z}; %#ok<NASGU>
  save2(sprintf([seg_dir, stack_config.image_prefix, ...
    evaluate_config.segmentation_suffix,'.mat'] , z), 'label_map');
end

fprintf('constructing linkage graph according to ground-truth assignment\n');
links_3D_ground_truth = {};
for z = 1:n_section-1
  label_pairs = unique([segmentation.wcat{z}(:), ...
    segmentation.wcat{z+1}(:)], 'rows');
  label_pairs = label_pairs(label_pairs(:,1)>0 & label_pairs(:,2)>0, :);
  links_3D_ground_truth{z} = [label_pairs, ...
    label_mapping_as_to_gt{z}(label_pairs(:,1)+1)== ...
    label_mapping_as_to_gt{z+1}(label_pairs(:,2)+1) & ...
    label_mapping_as_to_gt{z}(label_pairs(:,1)+1)~=0]; %#ok<AGROW>
end
  
config.linkage.train.manual_annotation_file = ...
  [image_dir, 'train_annotation_link_eval_k_fold.mat'];

% compute linkage confidences, evaluate them and also store in
% link_dump
% link_dump = struct('to_be_linked', struct([]), 'not_to_be_linked', struct([]));

n_threshold = length(evaluate_config.link_thresholds);

rand_score_false_split = zeros(n_threshold,1);
rand_score_false_split_bound = zeros(n_threshold,1);
rand_score_false_merge = zeros(n_threshold,1);
rand_score_false_merge_bound = zeros(n_threshold,1);

edit_score_n_split = zeros(n_threshold,1);
edit_score_n_merge = zeros(n_threshold,1);

test_brackets_min = [1, n_section-floor(n_section/2)+1];
test_brackets_max = [floor(n_section/2), n_section];

link_eval.links_3D_p = {};

for k = 1:2
  fprintf('-Fold %d/%d-\n', k, 2);
  test_bracket_min = test_brackets_min(k);
  test_bracket_max = test_brackets_max(k);
  
  %%%%%%%%%%%%%%%
  % Training
  %%%%%%%%%%%%%%%
  fprintf('Training ...\n');
  % retain only bodies belonging to training set.
  train_annotation = [];
  train_annotation.al = {};
  train_annotation.wcat = {};
  if(~isfield(evaluate_config, 'is_trained_on_manual_annotation') || ...
      ~evaluate_config.is_trained_on_manual_annotation)
    train_annotation.label_mapping_as_to_gt = {};
  end
  fprintf('Sections included for training:\n');
  for z = [1:(test_bracket_min-1), (test_bracket_max+1):n_section]
    fprintf('%d\n', z);
    train_annotation.al{end+1} = manual_annotation.al{z};
    if(isfield(evaluate_config, 'is_trained_on_manual_annotation') && ...
        evaluate_config.is_trained_on_manual_annotation)
      train_annotation.wcat{end+1} = manual_annotation.wcat{z};
    else
      train_annotation.wcat{end+1} = segmentation.wcat{z};
      train_annotation.label_mapping_as_to_gt{end+1} = ...
        label_mapping_as_to_gt{z};
    end
  end
  save2(config.linkage.train.manual_annotation_file, '-STRUCT', 'train_annotation');
  pipeline(710, config);
  fprintf('done.\n');
  
  %%%%%%%%%%%%%%%
  % Generate linkage graph for test sections
  %%%%%%%%%%%%%%%
  fprintf('Generating linkage graphs for test sections ...\n');
  config.stack.case_ids = test_bracket_min:test_bracket_max;
  pipeline(750, config);
  fprintf('done.\n');
  
  %%%%%%%%%%%%%%%
  % Evaluate linkage graph for test sections
  %%%%%%%%%%%%%%%
  fprintf('Evaluating for this fold ...\n');
  for z = test_bracket_min:(test_bracket_max-1)
    fprintf(sprintf('Loading linkage graph for section %d ...\n', z));
    linkage_graph = load2([get_reconstruction_dir(config), linkage_config.dir, ...
      'links_3D_p', '.', num2str(z), '.', model_config.type, '.', feature_config.type, ...
      config.linkage.apply.model_suffix, config.segmentation_choose.choice.seg_suffix, ...
      '.mat'], 'links_3D_p');
    links_3D_p = linkage_graph.links_3D_p{1,1};
    links_3D_p = links_3D_p(links_3D_p(:,1)>0 & links_3D_p(:,2)>0, :);
    link_eval.links_3D_p{z} = links_3D_p;
    
    % make segment ids non-overlaping
    segment_offset = max(find(label_mapping_as_to_gt{z})-1) + 1;
    links_3D_p(:,2) = links_3D_p(:,2) + segment_offset;
    
    lp = links_3D_ground_truth{z};
    lp(:,2) = lp(:,2) + segment_offset;
    sec_seg_2_body_map = get_sec_seg_2_body_map_from_links_3D_with_dummy(...
      {lp}, 0.5);
    seg_2_body_map_gt = sec_seg_2_body_map(:, 2:3);
%     seg_2_body_map_gt = ...
%       [ find(label_mapping_as_to_gt{z})-1, ...
%         label_mapping_as_to_gt{z}(label_mapping_as_to_gt{z}>0); ...
%         find(label_mapping_as_to_gt{z+1})-1 + segment_offset, ...
%         label_mapping_as_to_gt{z+1}(label_mapping_as_to_gt{z+1}>0)];
    
    % evaluate linkage confidences
    for t = 1:n_threshold
      linkage_threshold = evaluate_config.link_thresholds(t);
      fprintf('linkage_threshold: %g\n', linkage_threshold);
      
      % get the pixelwise body ids
      sec_seg_2_body_map = get_sec_seg_2_body_map_from_links_3D_with_dummy(...
        {links_3D_p}, linkage_threshold);
      seg_2_body_map = sec_seg_2_body_map(:,2:3);
      % find segments that are not in the linkage graph and each to a
      % distinct body
      absent_segment_id = seg_2_body_map_gt(...
        ~ismember(seg_2_body_map_gt(:,1), seg_2_body_map(:,1)), 1);
      absent_segment_body_offset = max(seg_2_body_map(:,2))+1;
      seg_2_body_map = [seg_2_body_map; absent_segment_id, ...
        (absent_segment_body_offset:absent_segment_body_offset+length(absent_segment_id)-1)']; %#ok<AGROW>
      
      fprintf('Computing rand score ...\n');
      
      rand_score = compute_rand_score(uint32(seg_2_body_map_gt), ...
        uint32(seg_2_body_map));
      rand_score_false_split(t) = rand_score_false_split(t) + rand_score(1);
      rand_score_false_split_bound(t) = rand_score_false_split_bound(t) + ...
        rand_score(2);
      rand_score_false_merge(t) = rand_score_false_merge(t) + rand_score(3);
      rand_score_false_merge_bound(t) = rand_score_false_merge_bound(t) + ...
        rand_score(4);
      
      fprintf('Estimating number of splits and mergers needed for correction ...\n');
      edit_score = compute_split_merge_error(seg_2_body_map_gt, seg_2_body_map);
      edit_score_n_split(t) = edit_score_n_split(t) + edit_score(1);
      edit_score_n_merge(t) = edit_score_n_merge(t) + edit_score(2);
      
      if(isfield(evaluate_config, 'display_links') && evaluate_config.display_links)
        %%%%%%%%%%%%%%%%
        % Displaying results visually
        %%%%%%%%%%%%%%%%
        
        H = fspecial('gaussian');
        bdry_labels = ones(size(manual_annotation.al{z}));
        bdry_labels(manual_annotation.wcat{z}==0) = 0;
        expanded_bdry_labels = ones(size(bdry_labels));
        expanded_bdry_labels(conv2(bdry_labels,H,'same')<1) = 0;
        centroids1 = regionprops(manual_annotation.wcat{z},'Centroid');
        first_image_wbdry = (255-double(manual_annotation.al{z}))...
          .*expanded_bdry_labels;
        first_image_wbdry(:,:,2) = 255-double(manual_annotation.al{z});
        first_image_wbdry(:,:,3) = (255-double(manual_annotation.al{z}))...
          .*expanded_bdry_labels;
        
        
        
        bdry_labels = ones(size(manual_annotation.al{z}));
        bdry_labels(manual_annotation.wcat{z+1}==0) = 0;
        expanded_bdry_labels = ones(size(bdry_labels));
        expanded_bdry_labels(conv2(bdry_labels,H,'same')<1) = 0;
        
        centroids2 = regionprops(manual_annotation.wcat{z+1},'Centroid');
        
        
        
        second_image_wbdry = (255-double(manual_annotation.al{z+1}))...
          .*expanded_bdry_labels;
        second_image_wbdry(:,:,2) = 255-double(manual_annotation.al{z+1});
        second_image_wbdry(:,:,3) = (255-double(manual_annotation.al{z+1}))...
          .*expanded_bdry_labels;
        
        offset = 255*ones(size(manual_annotation.al{z},1),20);
        
        Big_matrix = [first_image_wbdry(:,:,1),offset,second_image_wbdry(:,:,1)];
        Big_matrix(:,:,2) = [first_image_wbdry(:,:,2),offset,second_image_wbdry(:,:,2)];
        Big_matrix(:,:,3) = [first_image_wbdry(:,:,3),offset,second_image_wbdry(:,:,3)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure(1);
        axes('Position', [0 0 1 1]);
        set(gcf, 'Position', [20, 20, size(Big_matrix, 2), size(Big_matrix, 1)]);
        imshow(uint8(Big_matrix));title('Displaying True positive');
        
        true_positive = links_3D_p(links_3D_p(:,3)>= evaluate_config.link_thresholds(t)...
          &  links_3D_p(:,1)==links_3D_p(:,2),1:2);
        
        if size(true_positive,1)>=1
          
          junk = [centroids1(true_positive(1:end,1)).Centroid];
          
          if size(true_positive,1)==1
            st_pt = junk;
          else
            st_pt = (reshape(junk,[2,size(true_positive,1)]))';
          end
          
          junk = [centroids2(true_positive(1:end,2)).Centroid];
          
          if size(true_positive,1)==1
            end_pt = junk;
          else
            end_pt = (reshape(junk,[2,size(true_positive,1)]))';
          end
          
          end_pt(:,1) = end_pt(:,1)+size(manual_annotation.al{z},2)+20;
          line([(st_pt(:,1))';(end_pt(:,1))'],[(st_pt(:,2))';(end_pt(:,2))'],...
            'linewidth',3);
        else
          title('True positive Not Found');
        end
        
        saveas(1, [get_reconstruction_dir(config_0), config.linkage.dir, ...
          config.linkage.evaluate.save_dir, 'link_result.tp.', num2str(z), '.', ...
          model_config.type, '.', feature_config.type, ...
          config.linkage.apply.model_suffix, '.tif']);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure(2);imshow(uint8(Big_matrix));title('Displaying True Negative');
        
        true_negative = links_3D_p(links_3D_p(:,3)< evaluate_config.link_thresholds(t) & ...
          links_3D_p(:,1)~=links_3D_p(:,2),1:2);
        
        if size(true_negative,1)>=1
          
          junk = [centroids1(true_negative(1:end,1)).Centroid];
          
          if size(true_negative,1)==1
            st_pt = junk;
          else
            st_pt = (reshape(junk,[2,size(true_negative,1)]))';
          end
          
          junk = [centroids2(true_negative(1:end,2)).Centroid];
          
          if size(true_negative,1)==1
            end_pt = junk;
          else
            end_pt = (reshape(junk,[2,size(true_negative,1)]))';
          end
          
          end_pt(:,1) = end_pt(:,1)+size(manual_annotation.al{z},2)+20;
          line([(st_pt(:,1))';(end_pt(:,1))'],[(st_pt(:,2))';(end_pt(:,2))'],...
            'linewidth',3);
        else
          title('True Negative Not Found');
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(3); imshow(uint8(Big_matrix));title('Displaying False positive');
        
        false_positive = links_3D_p(links_3D_p(:,3)>= evaluate_config.link_thresholds(t) & ...
          links_3D_p(:,1)~=links_3D_p(:,2),1:2);
        
        if size(false_positive,1)>=1
          
          junk = [centroids1(false_positive(1:end,1)).Centroid];
          
          
          if size(false_positive,1)==1
            st_pt = junk;
          else
            st_pt = (reshape(junk,[2,size(false_positive,1)]))';
          end
          
          junk = [centroids2(false_positive(1:end,2)).Centroid];
          
          if size(false_positive,1)==1
            end_pt = junk;
          else
            end_pt = (reshape(junk,[2,size(false_positive,1)]))';
          end
          
          end_pt(:,1) = end_pt(:,1)+size(manual_annotation.al{z},2)+20;
          line([(st_pt(:,1))';(end_pt(:,1))'],[(st_pt(:,2))';(end_pt(:,2))'],...
            'linewidth',3);
        else
          title('False positive Not Found');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure(4); imshow(uint8(Big_matrix));title('Displaying False Negative');
        
        false_negative = links_3D_p(links_3D_p(:,3)< evaluate_config.link_thresholds(t) & ...
          links_3D_p(:,1)==links_3D_p(:,2),1:2);
        
        if size(false_negative,1)>=1
          
          junk = [centroids1(false_negative(1:end,1)).Centroid];
          
          if size(false_negative,1)==1
            st_pt = junk;
          else
            st_pt = (reshape(junk,[2,size(false_negative,1)]))';
          end
          
          junk = [centroids2(false_negative(1:end,2)).Centroid];
          
          if size(false_negative,1)==1
            end_pt = junk;
          else
            end_pt = (reshape(junk,[2,size(false_negative,1)]))';
          end
          end_pt(:,1) = end_pt(:,1)+size(manual_annotation.al{z},2)+20;
          line([(st_pt(:,1))';(end_pt(:,1))'],[(st_pt(:,2))';(end_pt(:,2))'],...
            'linewidth',3);
        else
          title('False negative Not Found');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end
    end
  end
  
end

link_eval.rand_score_false_split = rand_score_false_split;
link_eval.rand_score_false_split_bound = rand_score_false_split_bound;
link_eval.rand_score_false_merge = rand_score_false_merge;
link_eval.rand_score_false_merge_bound = rand_score_false_merge_bound;

link_eval.edit_score_n_split = edit_score_n_split;
link_eval.edit_score_n_merge = edit_score_n_merge;

figure;
hPlot = plot(rand_score_false_merge, rand_score_false_split, 's-');
axis([0, max(rand_score_false_merge_bound), ...
  0, max(rand_score_false_split_bound)]);
xlabel('false merge');
ylabel('false split');
hLegend = legend(hPlot, [model_config.type, ' ', feature_config.type, ' ', ...
  config.linkage.apply.model_suffix]);
set(hLegend, 'Interpreter', 'none');
title(['Linkage segment-wise rand ', model_config.type, ' ', ...
  feature_config.type, ' ', config.linkage.apply.model_suffix]);
check_for_dir([get_reconstruction_dir(config_0), config.linkage.dir, ...
  config.linkage.evaluate.save_dir]);
saveas(gcf, [get_reconstruction_dir(config_0), config.linkage.dir, ...
  config.linkage.evaluate.save_dir, 'link_seg_rand.', model_config.type, ...
  '.', feature_config.type, config.linkage.apply.model_suffix, '.fig']);

figure;
hPlot = plot(edit_score_n_split, edit_score_n_merge, 's-');
axis([0, max(max(edit_score_n_split), 1), 0, max(max(edit_score_n_merge), 1)]);
xlabel('number of required splits');
ylabel('number of required mergers');
hLegend = legend(hPlot, [model_config.type, ' ', feature_config.type, ' ', ...
  config.linkage.apply.model_suffix]);
set(hLegend, 'Interpreter', 'none');
title(['Linkage edit score ', model_config.type, ' ', ...
  feature_config.type, ' ', config.linkage.apply.model_suffix]);
check_for_dir([get_reconstruction_dir(config_0), config.linkage.dir, ...
  config.linkage.evaluate.save_dir]);
saveas(gcf, [get_reconstruction_dir(config_0), config.linkage.dir, ...
  config.linkage.evaluate.save_dir, 'link_edit.', model_config.type, ...
  '.', feature_config.type, config.linkage.apply.model_suffix, '.fig']);

save2([get_reconstruction_dir(config_0), config.linkage.dir, ...
  config.linkage.evaluate.save_dir, 'link_seg_rand_edit.', model_config.type, ...
  '.', feature_config.type, config.linkage.apply.model_suffix, '.mat'], 'link_eval');

% save2([get_reconstruction_dir(config), linkage_config.dir, ...
%   'link_dump_k_fold', '.', model_config.type, '.', feature_config.type, ...
%   config.linkage.apply.model_suffix, '.mat'], 'link_dump', '-v7.3');

return
end
